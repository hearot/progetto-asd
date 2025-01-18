#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <vector>

// #define INPUT_FILE "example.gfa"
// #define INPUT_FILE "example (cycle).gfa"
// #define INPUT_FILE "example (pdf).gfa"

// see https://github.com/pangenome/odgi/blob/master/test/DRB1-3123_unsorted.gfa
#define INPUT_FILE "DRB1-3123_unsorted.gfa"

// see https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms/chrY.hprc-v1.0-pggb.gfa.gz
// #define INPUT_FILE "chrY.pan.fa.a2fb268.4030258.6a1ecc2.smooth.gfa" // 197MB -- excluded in the repository

// see https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms/chrX.hprc-v1.0-pggb.gfa.gz
// #define INPUT_FILE "chrX.pan.fa.a2fb268.4030258.6a1ecc2.smooth.gfa" // 2.7GB -- excluded in the repository

// #define PATTERN "TTCA"
#define PATTERN "CTCTTGTAAGAAAAGTTCTCCAAGTCCCCACCCCACCCAGA"

#define SIGMA 131 // sigma for rolling hash (Karp-Rabin fingerprint)
#define PRIME 402309354485303 // random prime number for rolling hash

using namespace std;

class gfa_graph {
    protected:
        // from: [(to, (i_from, i_to)), ...]
        vector<vector<pair<int, pair<bool, bool>>>> adj;
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

        void add_edge(int from, bool i_from, int to, bool i_to) {
            adj[from].push_back({to, {i_from, i_to}});
            is_source[2 * to + i_to] = false;
        }

        gfa_graph get_acyclic() {
            gfa_graph G;

            remove_backward_edges_dfa(&G);
            return G;
        }

        string get_label(int i, bool sign = false) {
            return label[2*i + sign];
        }

        string get_label(pair<int, bool> node) {
            return label[2*node.first + node.second];
        }

        pair<int, bool> get_source() { // (i, sign)
            for (int i = 0; i < adj.size(); i++) {
                if (is_source[2 * i]) {
                    return {i, false};
                } else if (is_source[2 * i + 1]) {
                    return {i, true};
                }
            }

            return {-1, false}; // no source
        }

        pair<int, bool> get_dest(pair<int, bool> source_pair) { // (i, sign)
            bool is_dest = true;

            for (auto el : adj[source_pair.first]) {
                if (el.second.first == source_pair.second) {
                    is_dest = false;

                    auto node = get_dest({el.first, el.second.second});

                    if (node.first != -1) {
                        return node;
                    }
                }
            }

            if (is_dest) {
                return source_pair;
            } else {
                return {-1, false}; // no destination
            }
        }

        bool check_pattern(string P, pair<int, bool> source_pair, pair<int, bool> dest_pair) {
            string current = "";

            return traverse_and_find_pattern_if_acyclic(P, source_pair, dest_pair, &current);
        }
    private:
        bool contains_substring(string str, string substr) {
            int n = substr.size();
            int m = str.size();
            int H = 0;
            int Hp = 0; // pattern hash

            int SIGMAk = fast_exp_mod(SIGMA, n-1, PRIME); // SIGMA^(n-1) mod PRIME

            if (m < n) {
                return false;
            }

            int i = 0;

            for (; i < n; i++) {
                H = (SIGMA * H + str[i]) % PRIME;
                Hp = (SIGMA * Hp + substr[i]) % PRIME;
            }

            for(; i < m; i++) {
                if (H == Hp && str.compare(i - n, n, substr) == 0) { // same hash
                    return true;
                }

                H = (SIGMA * (H - SIGMAk * str[i - n]) + str[i]) % PRIME;
            }

            return false;
        }

        int fast_exp_mod(int b, int i, long long int p) {
            int result = 1;
            
            while (i) {
                if (i & 1) { // odd exponent check
                    result = (result * b) % p;
                }

                b = (b * b) % p;
                i >>= 1; // right shift to divide by 2
            }

            return result;
        }

        bool traverse_and_find_pattern_if_acyclic(string P, pair<int, bool> source_pair, pair<int, bool> dest_pair, string* current) { // no visited vectors since we assume the graph is acyclic            
            string current_label = get_label(source_pair);
            *current += current_label;

            if (source_pair == dest_pair) { // first then second coordinate eq check in lazy fashion
                if (contains_substring(*current, P)) { // true if the pattern is matched
                    return true;
                } else {
                    (*current).erase((*current).size() - current_label.size()); // remove label
                    return false;
                }
            } else {
                for (auto el : adj[source_pair.first]) {
                    if (el.second.first == source_pair.second) {
                        if (traverse_and_find_pattern_if_acyclic(P, {el.first, el.second.second}, dest_pair, current)) {
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
                visited[2*i] = false;
                visited[2*i+1] = false;
            }

            for (int i = 0; i < adj.size(); i++) {
                if (!visited[2*i]) {
                    traverse_remove_cycles(G, i, false); // traverse i+
                }

                if (!visited[2*i + 1]) {
                    traverse_remove_cycles(G, i, true); // traverse i-
                }
            }
        }

        void traverse_remove_cycles(gfa_graph* G, int s, bool sign) {
            in_stack[2*s + sign] = true;
            visited[2*s + sign] = true;

            for (auto p : adj[s]) {
                if (p.second.first != sign) {
                    continue;
                }

                if (!visited[2 * p.first + p.second.second]) {
                    G->add_edge(s, sign, p.first, p.second.second);
                    traverse_remove_cycles(G,  p.first, p.second.second);
                } else if (!in_stack[2 * p.first + p.second.second]) { // if it's not a backward edge
                    G->add_edge(s, sign, p.first, p.second.second);
                }
            }

            in_stack[2*s + sign] = false;
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
    ifstream input_file(INPUT_FILE);

    if (!input_file.is_open()) {
        return 1; // input file is already open somewhere
    }

    gfa_graph G;

    cout << "Reading input file \"" << INPUT_FILE << "\"..." << endl;

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

            G.add_edge(segments[tokens[0]], tokens[1][0] == '-', segments[tokens[2]], tokens[3][0] == '-');
        } while (getline(input_file, line));
    }

    input_file.close();
    cout << "Done!" << endl;

    cout << "Making the graph acyclic..." << endl;
    G = G.get_acyclic();
    cout << "Done!" << endl;

    cout << "Getting a source node..." << endl;
    auto source_pair = G.get_source();
    cout << "Done!" << endl;

    cout << "Getting a destination node..." << endl;
    auto dest_pair = G.get_dest(source_pair);
    cout << "Done!" << endl;

    cout << "Checking if pattern \"" << PATTERN "\" is contained within a path..." << endl;
    bool found_pattern = G.check_pattern(PATTERN, source_pair, dest_pair);

    if (found_pattern) {
        cout << "The pattern was found!" << endl;
    } else {
        cout << "The pattern was NOT found!" << endl;
    }

    return 0;
}
