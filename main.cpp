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

// see https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms/chrY.hprc-v1.0-pggb.gfa.gz
// #define DEFAULT_INPUT_FILE "chrY.pan.fa.a2fb268.4030258.6a1ecc2.smooth.gfa" // 197MB -- excluded in the repository

// see https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms/chrX.hprc-v1.0-pggb.gfa.gz
// #define DEFAULT_INPUT_FILE "chrX.pan.fa.a2fb268.4030258.6a1ecc2.smooth.gfa" // 2.7GB -- excluded in the repository

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

typedef unordered_map<int, int> hashmap;

struct gfa_node {
    int index;
    bool sign; // true (1) if negative, false (0) if positive
};

bool operator==(const gfa_node& first, const gfa_node& second) { // lazy
    if (first.index != second.index) {
        return false;
    }

    if (first.sign != second.sign) {
        return false;
    }

    return true;
}

class gfa_graph {
    protected:
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
        void add_segment(string segment) {
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

        void add_edge(gfa_node source, gfa_node dest) {
            adj[source.index].push_back({{dest.index, dest.sign}, source.sign});
            is_source[2 * dest.index + dest.sign] = false; // bool treated as int!
        }

        gfa_graph get_acyclic() {
            gfa_graph G;

            remove_backward_edges_dfa(&G);
            return G;
        }

        string get_label(int index, bool sign = false) {
            return label[2 * index + sign]; // bool treated as int!
        }

        string get_label(gfa_node node) {
            return label[2 * node.index + node.sign]; // bool treated as int!
        }

        gfa_node get_source() { // (i, sign)
            for (int i = 0; i < adj.size(); i++) {
                if (is_source[2 * i]) {
                    return {i, false};
                } else if (is_source[2 * i + 1]) {
                    return {i, true};
                }
            }

            return {-1, false}; // no source
        }

        gfa_node get_dest(gfa_node source) { // (i, sign)
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

        bool check_pattern(string pattern, gfa_node source, gfa_node dest) {
            string current = "";

            return traverse_and_find_pattern_if_acyclic(pattern, source, dest, &current);
        }

        void print_most_frequent_kmers(int len, int n, gfa_node source, gfa_node dest) { // (length of K-mers, top n)
            hashmap occ; // frequencies
            string current = "";

            traverse_and_count_occurrences_if_acyclic(source, dest, &current, occ, len);

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
        string convert_to_string(int u, int len) { // length is needed to distinguish between "ATA" and "TA" for example ('A' = 0)
            string res = "";
            int c = 0;
            int i = 0;

            while (u != 0) {
                i++;

                c = u % 4;
                u /= 4;

                switch (c) {
                    case 0:
                        res = "A" + res;
                        break;
                    case 1:
                        res = "T" + res;
                        break;
                    case 2:
                        res = "C" + res;
                        break;
                    case 3:
                        res = "G" + res;
                        break;
                }
            }

            if (i < len) {
                for (; i != len; i++) {
                    res = "A" + res;
                }
            }

            return res;
        }

        int convert_to_base_4(char c) {
            switch (c) {
                case 'A':
                    return 0;
                case 'T':
                    return 1;
                case 'C':
                    return 2;
                case 'G':
                    return 3;
                default:
                    return -1; // c not in {A, T, C, G}
            }
        }

        // similar to `contains_substring`, except it just counts occurrences
        // and doesn't do any pattern-search (i.e. doesn't actually use
        // the Karp-Rabin fingerprint)
        void count_occurrences(string str, hashmap& occ, int len) {
            int m = str.size();

            if (m < len) {
                return;
            }

            int SIGMAk = fast_exp(SIGMA, len - 1); // SIGMA^(len-1) with fast exp
            int s = 0;
            int i = 0;

            for (; i < len; i++) {
                s = s * 4 + convert_to_base_4(str[i]);
            }

            occ[s]++;

            for (; i < m; i++) {
                s = 4 * (s - convert_to_base_4(str[i-len]) * SIGMAk) + convert_to_base_4(str[i]);
                occ[s]++;
            }
        }

        bool contains_substring(string str, string substr) { // uses Karp-Rabin fingerprint
            int n = substr.size();
            int m = str.size();
            int H = 0;
            int Hp = 0; // pattern hash

            if (m < n) {
                return false;
            }

            int SIGMAk = fast_exp_mod(SIGMA, n-1, PRIME); // SIGMA^(n-1) mod PRIME with fast exp
            int i = 0;

            for (; i < n; i++) {
                H = (SIGMA * H + convert_to_base_4(str[i])) % PRIME;
                Hp = (SIGMA * Hp + convert_to_base_4(substr[i])) % PRIME;
            }

            for(; i < m; i++) {
                if (H == Hp && str.compare(i - n, n, substr) == 0) { // same hash
                    return true;
                }

                H = (SIGMA * (H - SIGMAk * convert_to_base_4(str[i - n])) + convert_to_base_4(str[i])) % PRIME;
            }

            return false;
        }

        int fast_exp(int base, int exp) {
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

        int fast_exp_mod(int base, int exp, long long int prime) {
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
        void traverse_and_count_occurrences_if_acyclic(gfa_node source, gfa_node final_dest, string* current, hashmap& occ, int len) { // no visited vectors since we assume the graph is acyclic
            string current_label = get_label(source);
            *current += current_label;

            if (source == final_dest) { // first then second coordinate eq check in lazy fashion
                count_occurrences(*current, occ, len);
            } else {
                for (auto edge : adj[source.index]) {
                    if (edge.source_sign == source.sign) {
                        traverse_and_count_occurrences_if_acyclic(edge.dest, final_dest, current, occ, len);
                    }
                }
            }

            (*current).erase((*current).size() - current_label.size()); // remove label
        }

        bool traverse_and_find_pattern_if_acyclic(string pattern, gfa_node source, gfa_node final_dest, string* current) { // no visited vectors since we assume the graph is acyclic            
            string current_label = get_label(source);
            *current += current_label;

            if (source == final_dest) { // first then second coordinate eq check in lazy fashion
                if (contains_substring(*current, pattern)) { // true if the pattern is matched
                    return true;
                } else {
                    (*current).erase((*current).size() - current_label.size()); // remove label
                    return false;
                }
            } else {
                for (auto edge : adj[source.index]) {
                    if (edge.source_sign == source.sign) {
                        if (traverse_and_find_pattern_if_acyclic(pattern, edge.dest, final_dest, current)) {
                            return true;
                        }
                    }
                }
            }

            (*current).erase((*current).size() - current_label.size()); // remove label
            return false;
        }

        void remove_backward_edges_dfa(gfa_graph* G) {
            for (int i = 0; i < adj.size(); i++) {
                G->add_segment(get_label(i));
                visited[2 * i] = false;
                visited[2 * i+1] = false;
            }

            for (int i = 0; i < adj.size(); i++) {
                if (!visited[2 * i]) {
                    traverse_remove_cycles(G, {i, false}); // traverse i+
                }

                if (!visited[2 * i + 1]) {
                    traverse_remove_cycles(G, {i, true}); // traverse i-
                }
            }
        }

        void traverse_remove_cycles(gfa_graph* G, gfa_node source) {
            in_stack[2 * source.index + source.sign] = true;
            visited[2 * source.index + source.sign] = true;

            for (auto edge : adj[source.index]) {
                if (edge.source_sign != source.sign) {
                    continue;
                }

                if (!visited[2 * edge.dest.index + edge.dest.sign]) {
                    G->add_edge(source, edge.dest);
                    traverse_remove_cycles(G, edge.dest);
                } else if (!in_stack[2 * edge.dest.index + edge.dest.sign]) { // if it's not a backward edge
                    G->add_edge(source, edge.dest);
                }
            }

            in_stack[2 * source.index + source.sign] = false;
        }

        string flip(string segment) {
            string flipped;
            flipped.reserve(segment.size()); // reserve space -- its length is the same as segment's

            for (char c : segment) {
                if (c == 'A') {
                    flipped = 'T' + flipped;
                } else if (c == 'T') {
                    flipped = 'A' + flipped;
                } else if (c == 'C') {
                    flipped = 'G' + flipped;
                } else {
                    flipped = 'C' + flipped;
                }
            }

            return flipped;
        }
};

int main() {
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

    cout << "Making the graph acyclic... ";
    G = G.get_acyclic();
    cout << "done!" << endl;

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
