#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

// #define DEFAULT_INPUT_FILE "example.gfa"
// #define DEFAULT_INPUT_FILE "example (cycle).gfa"
#define DEFAULT_INPUT_FILE "example (pdf).gfa"

// see https://github.com/pangenome/odgi/blob/master/test/DRB1-3123_unsorted.gfa
// #define DEFAULT_INPUT_FILE "DRB1-3123_unsorted.gfa"

// for input file "example (pdf).gfa"
#define DEFAULT_PATTERN "TTCA"

// for input file "DRB1-3123_unsorted.gfa"
// #define DEFAULT_PATTERN "CTCTTGTAAGAAAAGTTCTCCAAGTCCCCACCCCACCCAGA" 

// K-mer length
#define DEFAULT_K 3

// number of letters (A, T, C, G)
#define SIGMA 4 // sigma for rolling hash (Karp-Rabin fingerprint)
#define PRIME 402309354485303 // random prime number for rolling hash

using namespace std;

int BASE_4[44] = {0}; // 'A' = 41, 'T' = 84, 84-41=44 -- lookup table
char SIGMA_TABLE[4] = {'A', 'T', 'C', 'G'}; // reverse lookup table

typedef unordered_map<int, int> hashmap;

struct gfa_node {
    int index;
    bool sign; // true (1) if negative, false (0) if positive

    bool operator==(const gfa_node& other) const { // lazy
        return (index == other.index) && (sign == other.sign);
    }
};

class gfa_graph {
    private:
        struct edge {
            gfa_node dest;
            bool source_sign; // true if negative, zero if positive
        };

        // from: [{dest: {index, sign}, source_sign}, ...]
        vector<vector<edge>> adj;
        vector<string> label; // 2*i -> label of i+; 2*i+1 -> label of i-
        vector<bool> visited;
        vector<bool> in_stack;
        vector<bool> is_source;

    public:
        void add_segment(const string& segment) {
            adj.push_back({});
            label.push_back(segment);
            label.push_back(flip(segment));

            // i+
            visited.push_back(false);
            in_stack.push_back(false);
            is_source.push_back(true);

            // i-
            visited.push_back(false);
            in_stack.push_back(false);
            is_source.push_back(true);
        }

        void add_edge(const gfa_node& source, const gfa_node& dest) {
            adj[source.index].push_back({{dest.index, dest.sign}, source.sign});
            is_source[2 * dest.index + dest.sign] = false; // bool treated as int!
        }

        pair<gfa_graph, bool> get_acyclic() { // (graph, true) if the graph was cyclic, (graph, false) otherwhise 
            gfa_graph G;
            bool is_cyclic = false;

            remove_backward_edges_dfa(G, is_cyclic);
            return {G, is_cyclic};
        }

        string get_label(const gfa_node& node) const {
            return label[2 * node.index + node.sign]; // bool treated as int!
        }

        gfa_node get_source() const { // (i, sign)
            for (int i = 0; i < adj.size(); i++) {
                if (is_source[2 * i]) {
                    return {i, false};
                } else if (is_source[2 * i + 1]) {
                    return {i, true};
                }
            }

            return {-1, false}; // no source
        }

        gfa_node get_dest(const gfa_node& source) const { // (i, sign)
            bool is_dest = true;

            // edge -> (dest: (index, sign), source_sign)
            for (auto edge : adj[source.index]) {
                if (edge.source_sign == source.sign) {
                    is_dest = false;

                    gfa_node node = get_dest(edge.dest);

                    if (node.index != -1) {
                        return node;
                    }
                }
            }

            if (is_dest) {
                return source;
            } else {
                return {-1, false}; // no destination
            }
        }

        bool check_pattern(const string& pattern, const gfa_node& source, const gfa_node& dest) const {
            string current_str = "";
            int current_hash = 0;
            bool found = false;

            int pattern_hash = 0;
            int pattern_len = pattern.length();

            for (int i = 0; i < pattern_len; i++) { // compute pattern hash
                pattern_hash = (SIGMA * pattern_hash + BASE_4[pattern[i] - 'A']) % PRIME;
            }

            int SIGMAk = fast_exp_mod(SIGMA, pattern_len - 1, PRIME); // SIGMA^(pattern_len-1) mod PRIME with fast exp

            return traverse_and_find_pattern_if_acyclic(pattern, pattern.length(), pattern_hash, SIGMAk, source, dest, current_str, current_hash, found);
        }

        void print_most_frequent_kmers(const int& len, const int& n, gfa_node source, const gfa_node& dest) const { // (length of K-mers, top n)
            hashmap occ; // frequencies
            string current = "";

            traverse_and_count_occurrences_if_acyclic(source, dest, current, occ, len);

            // min-heap (el: (number of occ's, k-mer representation in base 4))
            priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> ranking;

            auto it = occ.begin();
            int i = 0;

            while (i < n && it != occ.end()) {
                ranking.push({it->second, it->first});

                it++;
                i++;
            }

            if (i == 0) {
                cout << "No k-mers were found!" << endl;
                return;
            }

            while (it != occ.end()) {
                if (it->second > ranking.top().first) {
                    ranking.pop();
                    ranking.push({it->second, it->first}); // min-heap of fixed size n
                }

                it++;
            }

            stack<pair<string, int>> top;
            pair<int, int> elem;

            while (!ranking.empty()) {
                elem = ranking.top();
                top.push({convert_to_string(elem.second, len), elem.first});
                ranking.pop();
            }

            i = 1;
            pair<string, int> top_elem;

            while (!top.empty()) {
                top_elem = top.top();
                cout << i << ". " << top_elem.first << " - " << top_elem.second << endl;
                top.pop();
                i++;
            }
        }

    private:
        string convert_to_string(int u, const int& len) const { // length is needed to distinguish between "ATA" and "TA" for example ('A' = 0)
            string res = string(len, 'A'); // pre-allocate with default value 'A' (= 0, in our base)

            int c;

            for (int i = 1; u != 0; i++) {
                c = u % 4;
                u /= 4;

                res[len - i] = SIGMA_TABLE[c];
            }

            return res;
        }

        // similar to `contains_substring`, except it just counts occurrences
        // and doesn't do any pattern-search (i.e. doesn't actually use
        // the Karp-Rabin fingerprint)
        void count_occurrences(string str, hashmap& occ, const int& len) const {
            int m = str.size();

            if (m < len) {
                return;
            }

            int SIGMAk = fast_exp(SIGMA, len - 1); // SIGMA^(len-1) with fast exp
            int s = 0;
            int i = 0;

            for (; i < len; i++) {
                s = s * 4 + BASE_4[str[i] - 'A'];
            }

            occ[s]++;

            for (; i < m; i++) {
                s = 4 * (s - BASE_4[str[i-len] - 'A'] * SIGMAk) + BASE_4[str[i] - 'A'];
                occ[s]++;
            }
        }

        int fast_exp(int base, int exp) const {
            int result = 1;
            
            while (exp) {
                if (exp & 1) { // odd exponent check
                    result = result * base;
                }

                base = base * base;
                exp >>= 1; // right shift to divide by 2
            }

            return result;
        }

        int fast_exp_mod(int base, int exp, const long long int& prime) const {
            int result = 1;
            
            while (exp) {
                if (exp & 1) { // odd exponent check
                    result = (result * base) % prime;
                }

                base = (base * base) % prime;
                exp >>= 1; // right shift to divide by 2
            }

            return result;
        }

        // similar to `traverse_and_find_pattern_if_acyclic`, but does not do pattern-searching;
        // instead it counts occurrences of K-mers.
        void traverse_and_count_occurrences_if_acyclic(gfa_node source, const gfa_node& final_dest, string& current, hashmap& occ, const int& len) const { // no visited vectors since we assume the graph is acyclic
            string current_label = get_label(source);
            current += current_label;

            if (source == final_dest) { // first then second coordinate eq check in lazy fashion
                count_occurrences(current, occ, len);
            } else {
                for (auto edge : adj[source.index]) {
                    if (edge.source_sign == source.sign) {
                        traverse_and_count_occurrences_if_acyclic(edge.dest, final_dest, current, occ, len);
                    }
                }
            }

            current.erase(current.size() - current_label.size()); // remove label
        }

        bool traverse_and_find_pattern_if_acyclic(const string& pattern, const int& pattern_len, const int& pattern_hash, const int& SIGMAk, gfa_node source, const gfa_node& final_dest, string current_str, int current_hash, bool found) const { // no visited vectors since we assume the graph is acyclic            
            if (!found) { // if the pattern hasn't been found yet, update the hash & the string to try matching it
                string current_label = get_label(source);

                int q = current_str.length();

                current_str += current_label;
                int n = current_str.length();

                int i = q;

                for (; i < pattern_len && i < n; i++) { // initial hash
                    current_hash = (SIGMA * current_hash + BASE_4[current_str[i] - 'A']) % PRIME;
                }

                if (n >= pattern_len) {
                    // first pattern matching
                    if (current_hash == pattern_hash && current_str.compare(i - pattern_len, pattern_len, pattern) == 0) { // same hash
                        found = true;
                    } else {
                        while (i < n) { // rolling hash with a sliding window
                            current_hash = (SIGMA * (current_hash - SIGMAk * BASE_4[current_str[i - pattern_len] - 'A']) + BASE_4[current_str[i] - 'A']) % PRIME;
                            i++;

                            if (current_hash == pattern_hash && current_str.compare(i - pattern_len, pattern_len, pattern) == 0) { // same hash
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }

            if (source == final_dest) { // first then second coordinate eq check in lazy fashion
                return found;
            } else {
                for (auto edge : adj[source.index]) {
                    if (edge.source_sign == source.sign) {
                        if (traverse_and_find_pattern_if_acyclic(pattern, pattern_len, pattern_hash, SIGMAk, edge.dest, final_dest, current_str, current_hash, found)) {
                            return true;
                        }
                    }
                }

                return false;
            }
        }

        void remove_backward_edges_dfa(gfa_graph& G, bool& is_cyclic) {
            for (int i = 0; i < adj.size(); i++) {
                G.add_segment(get_label({i, false}));
                visited[2 * i] = false;
                visited[2 * i+1] = false;
            }

            for (int i = 0; i < adj.size(); i++) {
                if (!visited[2 * i]) {
                    traverse_remove_cycles(G, {i, false}, is_cyclic); // traverse i+
                }

                if (!visited[2 * i + 1]) {
                    traverse_remove_cycles(G, {i, true}, is_cyclic); // traverse i-
                }
            }
        }

        void traverse_remove_cycles(gfa_graph& G, gfa_node source, bool& is_cyclic) {
            in_stack[2 * source.index + source.sign] = true;
            visited[2 * source.index + source.sign] = true;

            for (auto edge : adj[source.index]) {
                if (edge.source_sign != source.sign) {
                    continue;
                }

                if (!visited[2 * edge.dest.index + edge.dest.sign]) {
                    G.add_edge(source, edge.dest);
                    traverse_remove_cycles(G, edge.dest, is_cyclic);
                } else if (!in_stack[2 * edge.dest.index + edge.dest.sign]) { // if it's not a backward edge
                    G.add_edge(source, edge.dest);
                } else { // it's a backward edge, hence the original graph was cyclic
                    is_cyclic = true;
                }
            }

            in_stack[2 * source.index + source.sign] = false;
        }

        string flip(const string& segment) const {
            string flipped(segment.size(), '\0');

            int i = 0;
            auto c = segment.rbegin();

            while (c != segment.rend()) {
                if (*c == 'A') {
                    flipped[i] = 'T';
                } else if (*c == 'T') {
                    flipped[i] = 'A';
                } else if (*c == 'C') {
                    flipped[i] = 'G';
                } else {
                    flipped[i] = 'C';
                }

                c++;
                i++;
            }

            return flipped;
        }
};

int main() {
    BASE_4['A' - 'A'] = 0;
    BASE_4['T' - 'A'] = 1;
    BASE_4['C' - 'A'] = 2;
    BASE_4['G' - 'A'] = 3;

    string input_filename;

    cout << "Insert input file (default is \"" << DEFAULT_INPUT_FILE << "\"): ";
    getline(cin, input_filename);

    input_filename = input_filename == "" ? DEFAULT_INPUT_FILE : input_filename;

    ifstream input_file(input_filename);

    if (!input_file.is_open()) {
        cout << "Couldn't open the input file!" << endl;
        return 1; // couldn't open the input file
    }

    gfa_graph G;

    cout << endl << "Reading input file \"" << input_filename << "\"... ";

    {
        string line;
        string token;
        string segment_name;
        map<string, int> segments;
        int counter = 0;

        getline(input_file, line); // ignore the header

        while (getline(input_file, line)) { // segments
            if (line[0] == 'L') {
                break;
            }

            istringstream stream(line);

            getline(stream, token, '\t'); // S
            getline(stream, token, '\t');

            segment_name = token;

            getline(stream, token, '\t');

            G.add_segment(token);
            segments[segment_name] = counter;
            counter++;
        }

        do {
            if (line[0] != 'L') {
                break;
            }

            vector<string> tokens;
            istringstream stream(line);

            getline(stream, token, '\t'); // L

            for (int i = 0; i < 4; i++) {
                getline(stream, token, '\t');
                tokens.push_back(token);
            }

            G.add_edge({segments[tokens[0]], tokens[1][0] == '-'}, {segments[tokens[2]], tokens[3][0] == '-'});
        } while (getline(input_file, line));
    }

    input_file.close();
    cout << "done!" << endl;

    cout << "Checking whether the graph is acyclic... ";
    auto P = G.get_acyclic();
    cout << "done!" << endl;

    G = P.first; // G is now acyclic

    if (P.second) {
        cout << "The graph was cyclic, and has now been rendered acyclic!" << endl;
    } else {
        cout << "The graph was already acyclic!" << endl;
    }

    cout << "Getting a source node... ";
    auto source = G.get_source();
    cout << "done!" << endl;

    cout << "Getting a destination node... ";
    auto dest = G.get_dest(source);
    cout << "done!" << endl << endl;

    string pattern = "";

    cout << "Insert pattern to search for (default is \"" << DEFAULT_PATTERN << "\"): ";
    getline(cin, pattern);

    pattern = pattern == "" ? DEFAULT_PATTERN : pattern;

    cout << endl << "Checking if pattern \"" << pattern << "\" is contained within a path..." << endl;
    bool found_pattern = G.check_pattern(pattern, source, dest);

    if (found_pattern) {
        cout << "The pattern was found!" << endl << endl;
    } else {
        cout << "The pattern was NOT found!" << endl << endl;
    }

    int k;

    {
        string k_str;
        cout << "Insert a length for which the top 10 most frequent K-mers have to be found (default is " << DEFAULT_K << "): ";
        getline(cin, k_str);

        k = k_str == "" ? DEFAULT_K : stoi(k_str);
    }

    cout << endl << "Ranking the top 10 most frequent K-mers with K=" << k << ":" << endl;
    G.print_most_frequent_kmers(k, 10, source, dest);

    return 0;
}
