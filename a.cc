#include <cctype>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
using namespace std;

string ToUpper(string s) {
  string upper;
  for (char c : s) {
    upper += toupper(c);
  }

  return upper;
}

typedef unordered_set<string> StringSet;

bool ExistsInSet(const StringSet& set, const string& s) {
  auto it = set.find(s);
  return it != set.end();
}

void AddToSet(StringSet& set, const string& s) {
  assert(!ExistsInSet(set, s));
  set.insert(s);
}

struct Point {
  Point() {}
  Point(int r, int c) : row(r), col(c) {}

  friend ostream& operator<<(ostream& os, const Point& p);

  int row = 0;
  int col = 0;
};

struct Span {
  Span(Point p, int l, bool v) : point(p), len(l), vert(v) {}

  Point GetPoint(int i) const {
    assert(i >= 0 && i < len);
    if (vert) {
      return Point(point.row + i, point.col);
    } else {
      return Point(point.row, point.col + i);
    }
  }

  friend ostream& operator<<(ostream& os, const Span& s);

  Point point;
  int len;
  bool vert;
};

typedef vector<Span> Spans;

ostream& operator<<(ostream& os, const Span& s) {
  os << "(" << s.point << "," << s.len << ",";
  if (s.vert) {
    os << "vertical)";
  } else {
    os << "horizontal)";
  }
  return os;
}

ostream& operator<<(ostream& os, const Point& p) {
  os << "(" << p.row << "," << p.col << ")";
  return os;
}

struct Word {
  Word() {}
  Word(string s) : word(s) {}
  int len() const { return word.length(); }
  string word;
};

typedef vector<Word*> Words;
typedef unordered_map<string, Words> WordMap;

struct Attr {
  bool is_empty() const { return has_blanks && !has_letters; }
  bool is_partial() const { return has_blanks && has_letters; }
  bool is_full() const { return !has_blanks && has_letters; }
  bool has_letters = false;
  bool has_blanks = false;
};

struct Grid {
  Grid(string n) {
    name = n;
  }

  vector<string> lines;
  string name;
  Spans spans;

  int rows() const { return lines.size(); }
  int cols() const {
    if (lines.empty()) {
      return 0;
    }

    return lines[0].length();
  }

  int max_size() const { return max(rows(), cols()); }

  bool in_bounds(const Point& p) const {
    return (p.row < rows() && p.col < cols() && p.row >= 0 && p.col >= 0);
  }

  char box(const Point& p) const {
    assert(in_bounds(p));
    return lines[p.row][p.col];
  }

  bool is_block(const Point& p) const {
    return box(p) == '.';
  }

  bool is_blank(const Point& p) const {
    return box(p) == '-';
  }

  bool is_letter(const Point& p) const {
    char c = box(p);
    return c >= 'A' && c <= 'Z';
  }

  bool Next(Point& p, bool vert) const {
    if (vert) {
      p.row++;
      if (p.row >= rows()) {
        p.row = 0;
        p.col++;
      }
    } else {
      p.col++;
      if (p.col >= cols()) {
        p.col = 0;
        p.row++;
      }
    }
    return in_bounds(p);
  }

  string GetString(const Span& s, Attr& attr) const {
    string str;
    for (int i = 0; i < s.len; i++) {
      Point p = s.GetPoint(i);
      assert(in_bounds(p) && !is_block(p));
      if (is_letter(p)) {
        attr.has_letters = true;
      } else if (is_blank(p)) {
        attr.has_blanks = true;
      }
      str += box(p);
    }

    return str;
  }

  void write_box(const Point& p, char c) {
    assert(in_bounds(p));
    lines[p.row][p.col] = c;
  }

  void WriteString(const Span& s, const string& w) {
    assert(s.len == w.length());
    for (int i = 0; i < s.len; i++) {
      Point p = s.GetPoint(i);
      write_box(p, w[i]);
    }
  }

  void LoadFromFile(string filename) {
    ifstream f;
    f.open(filename);
    while (!f.eof()) {
      string line;
      getline(f, line);
      if (!line.empty() && line[0] != '#') {
        lines.push_back(line);
      }
    }
  }

  void Check() const {
    for (string s : lines) {
      assert(s.length() == cols());
    }
  }

  void Print() const {
    cout << "Grid: " << name << "\n"
         << "Row size is: " << rows() << "\n"
         << "Column size is: " << cols() << "\n"
         << "Max size is: " << max_size() << "\n";
    for (string s : lines) {
      cout << s << "\n";
    }
  }

  void FillSpans(bool vert) {
    Point p;
    Point start;
    bool span = false;
    int len = 0;
    bool eol;

    do {
      if (vert) {
        eol = (p.row + 1 == rows());
      } else {
        eol = (p.col + 1 == cols());
      }

      if (!is_block(p)) {
        if (!span) {
          span = true;
          start = p;
        }
        len++;
      }

      if (span && (eol || is_block(p))) {
        spans.push_back(Span(start, len, vert));
        span = false;
        len = 0;
      }
    } while(Next(p, vert));
  }

  void FillSpans() {
    assert(spans.empty());
    FillSpans(false);
    FillSpans(true);
  }

  void PrintSpans() const {
    for (const Span& s : spans) {
      Attr attr;
      cout << s << " " << GetString(s, attr) << "\n";
    }
  }
};

class Library {
  public:
    Library() {}
    ~Library() {
      for (Word* w : words_) {
        delete w;
      }
    }

    // Returns null if no matches found
    const Words* FindWord(const string& s) const {
      auto it = word_map_.find(s);
      if (it != word_map_.end()) {
        return &it->second;
      } else {
        return NULL;
      }
    }

    bool IsWord(string s) const {
      auto it = word_map_.find(s);
      if (it == word_map_.end()) {
        return false;
      }

      return true;
      // return word_map_.count(s) > 0;
    }

    void ComputeStats() {
      assert(counts_.empty());
      counts_.resize(18);
      for (const Word* w : words_) {
        int len = w->word.length();
        if (len < 18) {
          counts_[len]++;
        }
      }
    }

    void PrintStats() const {
      for (int i = 1; i < counts_.size(); i++) {
        cout << "length " << i << ": " << counts_[i] << "\n";
      }
    }

    string GetWord(int i) const {
      assert(i >= 0 && i < words_.size());
      return words_[i]->word;
    }

    void CreatePatternHash(Word* w) {
      int len = w->len();
      int num_patterns = 1 << len;
      for (int i = 0; i < num_patterns; i++) {
        string temp = w->word;
        for (int j = 0; j < len; j++) {
          if ((i >> j) & 1) {
            temp[j] = '-';
          }
        }
        word_map_[temp].push_back(w);
      }
    }

    void ReadFromFile(string filename, int max_size) {
      ifstream f;
      f.open(filename);
      while (f.is_open() && !f.eof()) {
        string line;
        getline(f, line);
        if (!line.empty()) {
          line = ToUpper(line);
          int len = line.length();
          if (line[len - 1] == '\r') {
            line = line.substr(0, len - 1);
          }
          if (line.length() <= max_size) {
            Word *w = new Word(line);
            words_.push_back(w);
            CreatePatternHash(w);
          }
        }
      }
    }

    void DebugBuckets() const {
      for (int i = 0; i < word_map_.bucket_count(); i++) {
        cout << i << ": " << word_map_.bucket_size(i) << "\n";
      }
    }
  private:
    Words words_;                // Master vector of words
    WordMap word_map_;           // Pattern hash ("D--" returns vector of <DAD, DIP, Dud, ...>)
    vector<int> counts_;
};

Library lib;

struct Slot {
  Slot(const Span& s, const string& p) : span(s), pattern(p) {}

  Span span;
  string pattern;
  friend ostream& operator<<(ostream& os, const Slot& s);
};

typedef vector<Slot> Slots;

ostream& operator<<(ostream& os, const Slot& s) {
  os << s.span << " " << s.pattern;
  return os;
}

class Solver {
public:
  Solver() {}
  void Solve(Grid& grid) {
    // cout << "Solving this grid:\n";
    Loop(grid, 0);
  }

private:
  void Loop(Grid grid, int depth) {
    depth++;
    // if (depth > 5) {
      // return;
    // }
    Slots full_slots;
    Slots partial_slots;
    Slots empty_slots;

    for (const Span &s : grid.spans) {
      Attr attr;
      string str = grid.GetString(s, attr);
      if (attr.is_empty()) {
        empty_slots.push_back(Slot(s, str));
      } else if (attr.is_full()) {
        full_slots.push_back(Slot(s, str));
      } else if (attr.is_partial()) {
        partial_slots.push_back(Slot(s, str));
      }
    }

    int num_empty = empty_slots.size();
    int num_full = full_slots.size();
    int num_partial = partial_slots.size();

    // cout << "empty = " << num_empty << "\n";
    // cout << "full = " << num_full << "\n";
    // cout << "partial = " << num_partial << "\n";

    for (const Slot& s : full_slots) {
      if (!lib.IsWord(s.pattern)) {
        // cout << s.pattern << " is not a word\n";
        return;
      }
    }

    StringSet set;
    for (const Slot& s : full_slots) {
      if (ExistsInSet(set, s.pattern)) {
        return;
      } else {
        AddToSet(set, s.pattern);
      }
    }

    if (num_partial == 0 && num_empty == 0) {
      cout << "Solution!!\n";
      grid.Print();
      return;
    }

    assert(num_partial > 0);
    CommitSlot(grid, partial_slots[0], depth);
  }

    void CommitSlot(Grid& grid, const Slot& slot, int depth) {
      // cout << "COMMIT slot " << slot << "\n";
      // cout << "Possible word choices: \n";
      const Words* words = lib.FindWord(slot.pattern);
      if (words) {
        for (const Word *w : *words) {
          grid.WriteString(slot.span, w->word);
          Loop(grid, depth);
        }
      } else {
        // cout << "No Matches\n";
      }
    }
};

int main() {
  Grid grid("My Grid");
  grid.LoadFromFile("test");
  grid.Check();
  grid.FillSpans();
  // grid.PrintSpans();

  lib.ReadFromFile("top_12000.txt", grid.max_size());

  Solver solver;
  solver.Solve(grid);
}
