#include <fstream>
#include <iostream>
#include <string>
#include <vector>


// Parameters
std::string IN;  // input benchmark file
std::string OUT; // output solution file
bool DEBUG;      // debug mode

// Common vars
int N; // grid size
int M; // number of pins


// Classes
struct Point;
using point_container_t = std::vector<Point>;
struct Edge;
using edge_container_t = std::vector<Edge>;
struct Cell;
using field_t = Cell;
struct Res;

// Functions
inline int wire_length(Point *a, Point *b);
inline int wire_length(int x1, int y1, int x2, int y2);
void read_params(int argc, char *argv[]);
void error(std::string &&msg);
point_container_t parse_bench(); // benchmark parser
int build_mst(point_container_t &pins, Res *res = nullptr);
int add_to_tree(Edge *edge, int new_id, int *tree, int m);
point_container_t hanan_grid(point_container_t &pins);
void route(field_t *field, point_container_t &pins);
void connect_pins(field_t *field, Point *p1, Point *p2);
inline void min_max(int &a, int &b, int &min, int &max);
void connect_x(field_t *field, int x, int y_min, int y_max);
void connect_y(field_t *field, int y, int x_min, int x_max);
void print_field(field_t *field);
int dump_solution(field_t *field, point_container_t &pins);


struct Point {
    int id;
    int x;
    int y;
    edge_container_t edges;

    Point() { }
    Point(int x_, int y_, int id_) : x(x_), y(y_), id(id_) { }
};

struct Edge {
    int id;
    Point *p1;
    Point *p2;
    int l;

    Edge() { }
    Edge(Point *p1_, Point *p2_, int id_, int l_) {
        p1 = p1_;
        p2 = p2_;
        id = id_;
        l = l_;
    }
    void init(Point *p1_, Point *p2_, int id_) {
        p1 = p1_;
        p2 = p2_;
        id = id_;
        l = wire_length(p1, p2);
    }
};

struct Cell {
    // pins
    static const char NO_PIN = '.';
    static const char PIN = '#';
    static const char HANAN = ',';
    static const char HANAN_PIN = '*';
    // wire
    static const char NO_WIRE = ' ';
    static const char M2 = '-';
    static const char M3 = '|';
    static const char M2M3 = '+';
    // wire halves
    static const uint8_t BEG = 0b01;
    static const uint8_t END = 0b10;
    static const uint8_t ALL = 0b11;

    uint8_t pin : 1;
    uint8_t hanan : 1;
    uint8_t m2 : 2; // horizontal (y = const)
    uint8_t m3 : 2; // vertical (x = const)

    Cell() : pin(0b0), hanan(0b0), m2(0b00), m3(0b00) { }
};

struct Res {
    int m;
    Edge *edges;
    int *tree;

    Res(int m_) : m(m_) {
        edges = new Edge[m * (m - 1) / 2]; // complete graph
        tree = new int[m];
    }
    ~Res() {
        delete[] edges;
        delete[] tree;
    }
};


int main(int argc, char *argv[]) {
    // read and init parameters
    read_params(argc, argv);

    // read pins from benchmark
    point_container_t pins = parse_bench();

    // build initial Minimum Spanning Tree, compute its wirelength
    int wl = build_mst(pins);

    // compute Hanan points
    point_container_t hanan = hanan_grid(pins);

    if (DEBUG) {
        std::cout << " Initial wirelength: " << wl << ", Hanan points: " << hanan.size() << '\n';
    }

    // Main Procedure: minimize wire length by adding Hanan points
    while (!hanan.empty()) {
        if (DEBUG) {
            std::cout << "  Hanan points: " << hanan.size() << '\n';
        }

        int hanan_points_checked = 0;
        int size = (int)pins.size();
        int x = -1, y = -1;
        pins.emplace_back(x, y, size); // reserve place for new Hanan point
        Res res(size + 1);
        // find one Hanan point which minimizes wire length the most
        for (auto &it : hanan) {
            // reset edges
            for (auto &iter : pins) {
                iter.edges.clear();
            }

            // add new Hanan point
            pins[size].x = it.x;
            pins[size].y = it.y;

            // build new tree
            int new_wl = build_mst(pins, &res);
            if (new_wl < wl) {
                wl = new_wl;
                x = it.x;
                y = it.y;
            }

            if (DEBUG) {
                if (!(++hanan_points_checked % 200)) {
                    std::cout << "   checked: " << hanan_points_checked << '/' << hanan.size() << '\n';
                }
            }
        }

        // exit if no better solutions
        if (x == -1) {
            pins.pop_back();
            break;
        }

        // build improved solution again
        for (auto &it : pins) {
            it.edges.clear();
        }
        pins[size].x = x;
        pins[size].y = y;
        ++size;
        build_mst(pins, &res);

        // remove useless Hanan points
        int checked = 0;
        for (int i = M; i < size - checked; ++i) {
            if (pins[i].edges.size() <= 2) {
                // try to remove this point
                int x_h = pins[i].x;
                int y_h = pins[i].y;
                int id_h = pins[i].id;
                int last = (int)pins.size() - 1;
                pins[i].x = pins[last].x;
                pins[i].y = pins[last].y;
                pins[i].id = pins[last].id;
                pins.pop_back();

                // check resulting tree
                for (auto &it : pins) {
                    it.edges.clear();
                }
                int new_wl = build_mst(pins);
                // it became worse, restore that point
                if (new_wl > wl) {
                    for (auto &it : pins) {
                        it.edges.clear();
                    }
                    pins.emplace_back(x_h, y_h, id_h);
                    build_mst(pins);
                    ++checked;
                    if (DEBUG) {
                        std::cout << " Removing Hanan point (" << x_h << ',' << y_h << ") results in increase of wire length: " << new_wl << " from " << wl << '\n';
                    }
                }
                else {
                    if (DEBUG) {
                        std::cout << " Removing Hanan point (" << x_h << ',' << y_h << "), new wire length: " << new_wl << " from " << wl << '\n';
                    }
                    wl = new_wl; // in case of we minimize
                }
            }
        }

        // remove added Hanan point
        auto hanan_end = hanan.end();
        for (auto it = hanan.begin(); it != hanan_end; ++it) {
            if ((it->x == x) && (it->y == y)) {
                hanan.erase(it);
                break;
            }
        }
    }
    if (DEBUG) {
        std::cout << " Resulting wire length: " << wl << '\n';
    }
    
    // create field
    field_t *field = new field_t[N * N];
    int m = (int)pins.size();
    for (int i = 0; i < m; ++i) {
        if (i >= M) {
            field[pins[i].x + N * pins[i].y].hanan = 1;
        }
        field[pins[i].x + N * pins[i].y].pin = 1;
    }

    // tracing
    for (auto &it : pins) {
        it.edges.clear();
    }
    build_mst(pins);
    route(field, pins);
    if (DEBUG) {
        print_field(field);
    }

    // dump solution to xml
    wl = dump_solution(field, pins);
    if (DEBUG) {
        std::cout << " Wire length after routing: " << wl << '\n';
    }

    delete[] field;
    return 0;
}


inline int wire_length(Point *a, Point *b) {
    return abs(a->x - b->x) + abs(a->y - b->y);
}
inline int wire_length(int x1, int y1, int x2, int y2) {
    return abs(x1 - x2) + abs(y1 - y2);
}

void read_params(int argc, char *argv[]) {
    // check number of arguments
    if ((argc != 3) && (argc != 4)) {
        error("Wrong amount of arguments: " + std::to_string(argc - 1));
    }

    // read required arguments
    IN = argv[1];

    OUT = argv[2];

    // read optional arguments
    if (argc == 4) {
        if ((argv[3][0] != 'D') || (argv[3][1] != '\0')) {
            error("Wrong argument: \"" + std::string(argv[3]) + "\". Did you mean \"D\"?");
        }
        else {
            DEBUG = true;
        }
    }
    else {
        DEBUG = false;
    }

    if (DEBUG) {
        std::cout << "Running with: IN=\"" << IN << "\", OUT=\"" << OUT << "\"; DEBUG mode is enabled.\n";
    }
}

void error(std::string &&msg) {
    std::cout << " !!! ERROR !!!\n"
              << msg << "\n\n"
              << "Use: \".\\VLSI_CAD_PA_STGEN IN OUT [D]\"\n"
              << " IN - input benchmark file\n"
              << " OUT - output solution file\n"
              << " D - (optional) debug mode\n\n";
    exit(EXIT_FAILURE);
}

int build_mst(point_container_t &pins, Res *res /*= nullptr*/) {
    // create edges
    int m = res ? res->m : (int)pins.size();
    int edges_num = m * (m - 1) / 2;
    Edge *edges = res ? res->edges : new Edge[edges_num];
    int cnt = 0;
    for (int i = 0; i < m - 1; ++i) {
        for (int j = i + 1; j < m; ++j) {
            edges[cnt].init(&pins[i], &pins[j], cnt++);
        }
    }

    // sort edges by length
    {
        auto cmp = [](void const *p1, void const *p2) { return ((Edge*)p1)->l - ((Edge*)p2)->l; };
        qsort(edges, edges_num, sizeof(Edge), cmp);
    }

    // build a tree
    int *tree = res ? res->tree : new int[m];
    for (int i = 0; i < m; ++i) {
        tree[i] = m + i;
    }
    cnt = 0;
    int wl = 0;
    for (int i = 0; i < edges_num; ++i) {
        int added_wl = add_to_tree(&edges[i], cnt, tree, m);
        if (added_wl != 0) {
            wl += added_wl;
            ++cnt;
            if (cnt == m - 1) {
                break;
            }
        }
    }

    if (!res) {
        delete[] tree;
        delete[] edges;
    }
    return wl;
}

int add_to_tree(Edge *edge, int new_id, int *tree, int m) {
    int id1 = tree[edge->p1->id];
    int id2 = tree[edge->p2->id];

    // cycle
    if (id1 == id2) {
        return 0;
    }

    if (id1 >= m) {
        if (id2 >= m) { // p1 new, p2 new
            tree[edge->p1->id] = new_id;
            tree[edge->p2->id] = new_id;
        }
        else { // p1 new
            tree[edge->p1->id] = id2;
        }
    }
    else {
        if (id2 >= m) { // p2 new
            tree[edge->p2->id] = id1;
        }
        else {
            for (int i = 0; i < m; ++i) {
                if (tree[i] == id2) {
                    tree[i] = id1;
                }
            }
        }
    }

    edge->p1->edges.emplace_back(edge->p1, edge->p2, new_id, edge->l);
    edge->p2->edges.emplace_back(edge->p1, edge->p2, new_id, edge->l);
    return edge->l;
}

point_container_t hanan_grid(point_container_t &pins) {
    // collect Xs and Ys with pins
    bool *xs = new bool[N];
    for (int i = 0; i < N; ++i) {
        xs[i] = false;
    }
    bool *ys = new bool[N];
    for (int i = 0; i < N; ++i) {
        ys[i] = false;
    }
    for (auto &it : pins) {
        xs[it.x] = true;
        ys[it.y] = true;
    }

    point_container_t ret;
    for (int y = 0; y < N; ++y) {
        if (ys[y]) {
            int yN = y * N;
            for (int x = 0; x < N; ++x) {
                if (xs[x]) {
                    bool exist = false;
                    for (auto &it : pins) {
                        if ((x == it.x) && (y == it.y)) {
                            exist = true;
                            break;
                        }
                    }
                    if (!exist) {
                        ret.emplace_back(x, y, -1); // add point w/o ID
                    }
                }
            }
        }
    }

    delete[] ys;
    delete[] xs;
    return ret;
}

void route(field_t *field, point_container_t &pins) {
    // Depth-First Search
    int m = (int)pins.size();
    Point **stack = new Point*[m];
    bool *e_visited = new bool[m - 1];
    for (int i = 0; i < m - 1; ++i) {
        e_visited[i] = false;
    }

    // push
    stack[0] = &(pins[0]);
    int elem = 1;

    while (elem != 0) {
        // pop
        Point *p = stack[elem - 1];
        --elem;

        if (DEBUG) {
            if (p->edges.size() > 4) {
                std::cout << " Point x=" << p->x << " y=" << p->y << " has more than 4 edges (" << p->edges.size() << ")\n";
            }
        }

        for (auto &it : p->edges) {
            if (!e_visited[it.id]) {
                e_visited[it.id] = true;
                Point *another = (it.p1 != p) ? it.p1 : it.p2;

                // connect pins
                connect_pins(field, p, another);

                // push
                stack[elem] = another;
                ++elem;
            }
        }
    }

    delete[] e_visited;
    delete[] stack;
}

void connect_pins(field_t *field, Point *p1, Point *p2) {
    // constant X: m3 level
    if (p1->x == p2->x) {
        int y_min, y_max;
        min_max(p1->y, p2->y, y_min, y_max);
        connect_x(field, p1->x, y_min, y_max);
    }
    // constant Y: m2 level
    else if (p1->y == p2->y) {
        int x_min, x_max;
        min_max(p1->x, p2->x, x_min, x_max);
        connect_y(field, p1->y, x_min, x_max);
    }
    // X & Y
    else {
        // another verticles of the rectangle
        int x1 = p1->x;
        int y1 = p2->y;
        int x2 = p2->x;
        int y2 = p1->y;

        // chose the closest to a corner
        int x = x1;  // better
        int y = y1;
        int x0 = x2; // connect to pin
        int y0 = y2;
        {
            int wl = wire_length(x1, y1, 0, 0);
            int wl_ = wire_length(x1, y1, 0, N - 1);
            if (wl_ < wl) {
                wl = wl_;
            }
            wl_ = wire_length(x1, y1, N - 1, 0);
            if (wl_ < wl) {
                wl = wl_;
            }
            wl_ = wire_length(x1, y1, N - 1, N - 1);
            if (wl_ < wl) {
                wl = wl_;
            }

            bool remap = false;
            wl_ = wire_length(x2, y2, 0, 0);
            if (wl_ < wl) {
                wl = wl_;
                remap = true;
            }
            wl_ = wire_length(x2, y2, 0, N - 1);
            if (wl_ < wl) {
                wl = wl_;
                remap = true;
            }
            wl_ = wire_length(x2, y2, N - 1, 0);
            if (wl_ < wl) {
                wl = wl_;
                remap = true;
            }
            wl_ = wire_length(x2, y2, N - 1, N - 1);
            if (wl_ < wl) {
                wl = wl_;
                remap = true;
            }
            if (remap) {
                x = x2;
                y = y2;
                x0 = x1;
                y0 = y1;
            }
        }

        // connect via better way
        int x_min, x_max, y_min, y_max;
        min_max(y0, y, y_min, y_max);
        min_max(x0, x, x_min, x_max);
        connect_x(field, x, y_min, y_max);
        connect_y(field, y, x_min, x_max);
    }
}

inline void min_max(int &a, int &b, int &min, int &max) {
    if (a < b) {
        min = a;
        max = b;
    }
    else {
        min = b;
        max = a;
    }
}

void connect_x(field_t *field, int x, int y_min, int y_max) {
    field[x + y_min * N].m3 |= Cell::BEG;
    field[x + y_max * N].m3 |= Cell::END;
    for (int y = y_min + 1; y <= y_max - 1; ++y) {
        if (DEBUG) {
            if (field[x + y * N].m3) {
                std::cout << "Connect X: x=" << x << " y=[" << y_min << ';' << y_max << "], wire exists at M3 at y=" << y << '\n';
            }
            if (field[x + y * N].pin) {
                std::cout << "Connect X: x=" << x << " y=[" << y_min << ';' << y_max << "], met pin at y=" << y << '\n';
            }
            if (field[x + y * N].hanan) {
                std::cout << "Connect X: x=" << x << " y=[" << y_min << ';' << y_max << "], met hanan point at y=" << y << '\n';
            }
        }
        field[x + y * N].m3 = Cell::ALL;
    }
}

void connect_y(field_t *field, int y, int x_min, int x_max) {
    int yN = y * N;
    field[x_min + yN].m2 |= Cell::BEG;
    field[x_max + yN].m2 |= Cell::END;
    for (int x = x_min + 1; x <= x_max - 1; ++x) {
        if (DEBUG) {
            if (field[x + yN].m2) {
                std::cout << "Connect Y: x=[" << x_min << ';' << x_max << "] y=" << y << ", wire exists at M2 at x=" << x << '\n';
            }
            if (field[x + yN].pin) {
                std::cout << "Connect Y: x=[" << x_min << ';' << x_max << "] y=" << y << ", met pin at x=" << x << '\n';
            }
            if (field[x + yN].hanan) {
                std::cout << "Connect Y: x=[" << x_min << ';' << x_max << "] y=" << y << ", met hanan point at x=" << x << '\n';
            }
        }
        field[x + yN].m2 = Cell::ALL;
    }
}

void print_field(field_t *field) {
    for (int y = 0; y < N; ++y) {
        int yN = y * N;

        // print pins
        for (int x = 0; x < N; ++x) {
            char c = 0;
            if (field[x + yN].pin && field[x + yN].hanan) {
                c = Cell::HANAN_PIN;
            }
            else if (field[x + yN].pin) {
                c = Cell::PIN;
            }
            else if (field[x + yN].hanan) {
                c = Cell::HANAN;
            }
            else {
                c = Cell::NO_PIN;
            }
            std::cout << ' ' << c;
        }

        std::cout << "    ";
        // print pins|wire
        for (int x = 0; x < N; ++x) {
            char c = 0;
            // pins
            if (field[x + yN].pin && field[x + yN].hanan) {
                c = Cell::HANAN_PIN;
            }
            else if (field[x + yN].pin) {
                c = Cell::PIN;
            }
            else if (field[x + yN].hanan) {
                c = Cell::HANAN;
            }
            else {
                c = Cell::NO_PIN;
            }
            std::cout << c;

            c = 0;
            // wires
            if (field[x + yN].m2 && field[x + yN].m3) {
                c = Cell::M2M3;
            }
            else if (field[x + yN].m2) {
                c = Cell::M2;
            }
            else if (field[x + yN].m3) {
                c = Cell::M3;
            }
            else {
                c = Cell::NO_WIRE;
            }
            std::cout << c;
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}

int dump_solution(field_t *field, point_container_t &pins) {
    std::ofstream solution(OUT);
    if (solution.is_open()) {
        // header
        // "<net grid_size="50" pin_count="30">"
        solution << "<net grid_size=\"" << N << "\" pin_count=\"" << M << "\">\n";

        for (int i = 0; i < M; ++i) {
            // input pins
            // "    <point x="3" y="3" layer="pins" type="pin" />"
            solution << "    <point x=\"" << pins[i].x << "\" y=\"" << pins[i].y << "\" layer=\"pins\" type=\"pin\" />\n";
            // via from pins to m2 (y = const, horizontal)
            // "    <point x="3" y="3" layer="pins_m2" type="via" />"
            solution << "    <point x=\"" << pins[i].x << "\" y=\"" << pins[i].y << "\" layer=\"pins_m2\" type=\"via\" />\n";
        }

        int NN = N * N;
        for (int i = 0; i < NN; ++i) {
            if (field[i].m3 && (field[i].m2 || field[i].pin)) {
                // via from m2 (y = const, horizontal) to m3 (x = const, vertical)
                // "    <point x="3" y="3" layer="m2_m3" type="via" />"
                solution << "    <point x=\"" << (i % N) << "\" y=\"" << (i / N) << "\" layer=\"m2_m3\" type=\"via\" />\n";
            }
        }

        int wl = 0;
        for (int y = 0; y < N; ++y) {
            int yN = y * N;
            for (int x = 0; x < N; ++x) {
                if (field[x + yN].m2 == Cell::BEG) {
                    int y1 = y, y2 = y;
                    int x1 = x, x2 = x;
                    while (x2 != N - 1) {
                        if (field[x2 + 1 + yN].m2 == Cell::ALL) {
                            ++x2;
                        }
                        else if (field[x2 + 1 + yN].m2 == Cell::END) {
                            ++x2;
                            break;
                        }
                        else {
                            if (DEBUG) {
                                std::cout << " Y-segment with unexpected end at: x=" << (x2 + 1) << " y=" << y << '\n';
                            }
                            break;
                        }
                    }
                    if (x2 != x1) {
                        // Y-segment (horizontal) on m2 (different X)
                        // "    <segment x1="1" y1="2" x2="3" y2="2" layer="m2" />
                        solution << "    <segment x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" layer=\"m2\" />\n";
                        wl += (x2 - x1);
                        x = x2;
                    }
                    else {
                        if (DEBUG) {
                            std::cout << " Y-segment with 0 length: x=" << x << " y=" << y << '\n';
                        }
                    }
                }
            }
        }
        for (int x = 0; x < N; ++x) {
            for (int y = 0; y < N; ++y) {
                if (field[x + y * N].m3 == Cell::BEG) {
                    int x1 = x, x2 = x;
                    int y1 = y, y2 = y;
                    while (y2 != N - 1) {
                        if (field[x + (y2 + 1) * N].m3 == Cell::ALL) {
                            ++y2;
                        }
                        else if (field[x + (y2 + 1) * N].m3 == Cell::END) {
                            ++y2;
                            break;
                        }
                        else {
                            if (DEBUG) {
                                std::cout << " X-segment with unexpected end at: x=" << x << " y=" << (y2 + 1) << '\n';
                            }
                            break;
                        }
                    }
                    if (y2 != y1) {
                        // X-segment (vertical) on m3 (different Y)
                        // "    <segment x1="2" y1="1" x2="2" y2="3" layer="m3" />
                        solution << "    <segment x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" layer=\"m3\" />\n";
                        wl += (y2 - y1);
                        y = y2;
                    }
                    else {
                        if (DEBUG) {
                            std::cout << " X-segment with 0 length: x=" << x << " y=" << y << '\n';
                        }
                    }
                }
            }
        }

        // length
        // "<length v="1" />"
        solution << "<length v=\"" << wl << "\" />\n";

        // footer
        // "</net>"
        solution << "</net>\n";

        solution.close();
        return wl;
    }
    else {
        error("Couldn't open solution file \"" + OUT + '\"');
    }
    return 0;
}


// -------- Benchmark file parser --------
const int MAX_ATTR_NAME_LEN = 10;

struct Attr;

void parse_header(std::ifstream &in);
void read_word(std::ifstream &in, const char *word);
void error_parse(std::string &&msg);
void read_attr(std::ifstream &in, Attr *attr, int attr_num);
void read_word_till(std::ifstream &in, char *buf, char till);
void parse_pin(std::ifstream &in, int *x, int *y);
void parse_footer(std::ifstream &in);

struct Attr {
    const char *name;
    int *parse_value;
    const char *ignore_value;
    bool seen;

    Attr(const char *name_, int *parse) {
        name = name_;
        parse_value = parse;
        //ignore_value = nullptr; // will not be accessed
        seen = false;
    }
    Attr(const char *name_, const char *ignore) {
        name = name_;
        parse_value = nullptr;
        ignore_value = ignore;
        seen = false;
    }
};

point_container_t parse_bench() {
    std::ifstream bench_file(IN);
    if (bench_file.is_open()) {
        // <net grid_size="50" pin_count="30">
        parse_header(bench_file);

        point_container_t ret;
        // <point x="3" y="3" layer="pins" type="pin" />
        for (int i = 0; i < M; ++i) {
            int x, y;
            parse_pin(bench_file, &x, &y);
            ret.emplace_back(x, y, i);
        }

        // </net>
        parse_footer(bench_file);

        bench_file.close();
        return ret;
    }
    else {
        error("Couldn't open benchmark file \"" + IN + '\"');
    }
    return point_container_t();
}

void parse_header(std::ifstream &in) {
    read_word(in, "<net");

    Attr attr[2] = { { "grid_size", &N }, { "pin_count", &M } };
    read_attr(in, attr, 2);

    read_word(in, ">");
}

void read_word(std::ifstream &in, const char *word) {
    while (true) {
        int c = in.get();
        if (c == EOF) {
            error_parse("Unexpected EOF");
        }
        if (isspace(c)) {
            continue;
        }
        if (*word != c) {
            error_parse("Unexpected beginning: \"" + std::string(1, (char)c) + "\" of the word \"" + std::string(word) + '\"');
        }
        ++word;
        break;
    }
    while (*word) {
        int c = in.get();
        if (c == EOF) {
            error_parse("Unexpected EOF");
        }
        if (*word != c) {
            error_parse("Unexpected char: \"" + std::string(1, (char)c) + "\" instead of \"" + std::string(word) + '\"');
        }
        ++word;
    }
}

void error_parse(std::string &&msg) {
    std::cout << " !!! BENCHMARK PARSING ERROR !!!\n"
              << msg << "\n\n";
    exit(EXIT_FAILURE);
}

void read_attr(std::ifstream &in, Attr *attr, int attr_num) {
    char buf[MAX_ATTR_NAME_LEN];
    for (int i = 0; i < attr_num; ++i) {
        read_word_till(in, buf, '=');

        read_word(in, "=");
        read_word(in, "\"");

        bool found = false;
        for (int j = 0; j < attr_num; ++j) {
            if (!strcmp(attr[j].name, buf)) {
                if (attr[j].seen) {
                    error_parse("Multiple attribute definition: \"" + std::string(attr[j].name) + '\"');
                }
                attr[j].seen = true;
                if (attr[j].parse_value) {
                    in >> *(attr[j].parse_value);
                }
                else {
                    read_word(in, attr[j].ignore_value);
                }
                found = true;
                break;
            }
        }
        if (!found) {
            error_parse("Unknown attribute: \"" + std::string(buf) + '\"');
        }

        read_word(in, "\"");
    }
}

void read_word_till(std::ifstream &in, char *buf, char till) {
    while (true) {
        int c = in.get();
        if (c == EOF) {
            error_parse("Unexpected EOF");
        }
        if (isspace(c)) {
            continue;
        }
        buf[0] = (char)c;
        break;
    }
    int pos = 1;
    while (true) {
        int c = in.get();
        if (c == EOF) {
            error_parse("Unexpected EOF");
        }
        if (((char)c == till) || isspace(c) || (pos == MAX_ATTR_NAME_LEN - 1)) {
            in.putback((char)c);
            buf[pos] = '\0';
            break;
        }
        buf[pos] = (char)c;
        ++pos;
    }
}

void parse_pin(std::ifstream &in, int *x, int *y) {
    read_word(in, "<point");

    Attr attr[4] = { { "x", x }, { "y", y }, { "layer", "pins" }, { "type", "pin" } };
    read_attr(in, attr, 4);

    read_word(in, "/>");
}

void parse_footer(std::ifstream &in) {
    read_word(in, "</net");

    read_word(in, ">");
}

