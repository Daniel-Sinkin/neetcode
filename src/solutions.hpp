// src/solutions.hpp
#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <print>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ds_neetcode {

using usize = std::size_t;

using lc_chars = std::array<int, 26>; // lower case ascii chars

struct TreeNode {
    int val;
    TreeNode *left;
    TreeNode *right;
    TreeNode() : val(0), left(nullptr), right(nullptr) {}
    TreeNode(int x) : val(x), left(nullptr), right(nullptr) {}
    TreeNode(int x, TreeNode *left, TreeNode *right) : val(x), left(left), right(right) {}
};

constexpr bool test_problem_1 = false;
constexpr bool test_problem_2 = false;
constexpr bool test_problem_3 = false;
constexpr bool test_problem_4 = false;
constexpr bool test_problem_5 = false;
constexpr bool test_problem_6 = true;

bool contains_duplicates(const std::vector<int> &nums) {
    std::unordered_set<int> my_set{};
    my_set.reserve(nums.size());

    for (const int elem : nums) {
        if (my_set.contains(elem)) {
            return true;
        }
        my_set.insert(elem);
    }
    return false;
}

void problem_1_test() {
    using Input = std::vector<int>;
    using Output = bool;
    using Case = std::pair<Output, Input>;
    std::vector<Case> test_cases{};
    test_cases.push_back({true, {1, 2, 3, 3}});
    test_cases.push_back({false, {1, 2, 3, 4}});

    for (const auto &[target, input] : test_cases) {
        assert(contains_duplicates(input) == target);
    }
}

bool is_anagram(const std::string &s, const std::string &t) {
    if (s.size() != t.size()) {
        return false;
    }

    lc_chars seen_chars{};

    for (const char sc : s) {
        seen_chars[static_cast<usize>(sc - 'a')] += 1;
    }

    for (const char tc : t) {
        seen_chars[static_cast<usize>(tc - 'a')] -= 1;
    }

    for (usize i = 0; i < 26; ++i) {
        if (seen_chars[i] != 0) {
            return false;
        }
    }

    return true;
}

void problem_2_test() {
    using Input = std::pair<std::string, std::string>;
    using Output = bool;
    using Case = std::pair<Output, Input>;
    std::vector<Case> test_cases{};
    test_cases.push_back({true, {"racecar", "carrace"}});
    test_cases.push_back({false, {"jar", "jam"}});

    for (const auto &[target, input] : test_cases) {
        assert(is_anagram(input.first, input.second) == target);
    }
}

std::vector<int> two_sum(const std::vector<int> &nums, int target) {
    std::unordered_map<int, usize> my_map{};
    for (usize i{0zu}; i < nums.size(); ++i) {
        const int num{nums[i]};
        auto it = my_map.find(target - num);
        if (it != my_map.end()) {
            const auto &[k, idx_found] = *it;
            (void)k;
            const int idx0 = static_cast<int>(i);
            const int idx1 = static_cast<int>(idx_found);
            return {std::min(idx0, idx1), std::max(idx0, idx1)};
        }
        my_map.insert({num, i});
    }
    std::unreachable();
}

void problem_3_test() {
    using Input = std::pair<std::vector<int>, int>;
    using Output = std::vector<int>;
    using Case = std::pair<Output, Input>;
    std::vector<Case> test_cases{};
    test_cases.push_back({{0, 1}, {{3, 4, 5, 6}, 7}});
    test_cases.push_back({{0, 2}, {{4, 5, 6}, 10}});
    test_cases.push_back({{0, 1}, {{5, 5}, 10}});

    for (const auto &[target, input] : test_cases) {
        assert(two_sum(input.first, input.second) == target);
    }
}

struct LCCharsHash {
    usize operator()(const lc_chars &a) const noexcept {
        usize h{0zu};
        for (usize i{0zu}; i < 26; ++i) {
            h = 31 * h + static_cast<usize>(a[i]);
        }
        return h;
    }
};

static lc_chars chars_counter(const std::string &s) {
    lc_chars retval{};
    for (const char c : s) {
        retval[static_cast<usize>(c - 'a')] += 1;
    }
    return retval;
}

std::vector<std::vector<std::string>>
group_anagrams(const std::vector<std::string> &strs) {
    std::unordered_map<lc_chars, std::vector<std::string>, LCCharsHash> groups{};
    groups.reserve(strs.size());

    for (const auto &s : strs) {
        groups[chars_counter(s)].push_back(s);
    }

    std::vector<std::vector<std::string>> out{};
    out.reserve(groups.size());
    for (auto &[k, v] : groups) {
        (void)k;
        out.push_back(std::move(v));
    }
    return out;
}

void problem_4_test() {
    using Input = std::vector<std::string>;
    using Output = std::vector<std::vector<std::string>>;
    using Case = std::pair<Output, Input>;

    std::vector<Case> test_cases{};

    test_cases.push_back({
        {{"hat"}, {"pots", "tops", "stop"}, {"act", "cat"}},
        {"act", "pots", "tops", "cat", "stop", "hat"},
    });

    test_cases.push_back({
        {{"x"}},
        {"x"},
    });

    test_cases.push_back({
        {{""}},
        {""},
    });

    for (auto &[expected, input] : test_cases) {
        assert(group_anagrams(input) == expected);
    }
}

void problem_5_test() {
}

namespace Problem6 {
using Decoded = std::vector<std::string>;
using Encoded = std::string;

Encoded encode(const Decoded &strs) {
    if (strs.empty()) {
        return "";
    }
    std::string out;
    for (const std::string &word : strs) {
        const std::string word_size_s = std::to_string(word.size());
        if (word_size_s.size() == 1) {
            out.push_back('1');
        } else if (word_size_s.size() == 2) {
            out.push_back('2');
        } else {
            out.push_back('3');
        }
        out.append(word_size_s);
        out.append(word);
    }
    return out;
}

Decoded decode(const Encoded &s) {
    Decoded out;
    usize i{0zu};
    while (i < s.size()) {
        const usize read_digit = static_cast<usize>(s[i] - '0');
        ++i;
        const usize n_chars = static_cast<usize>(std::stoul(s.substr(i, read_digit)));
        i += read_digit;
        out.push_back(std::string{s.substr(i, n_chars)});
        i += n_chars;
    }
    return out;
}

void problem_6_test() {
    using Input = Decoded;
    const Input input = {"1Hello,", "my", "old", "Friend!", ""};
    std::string encoded = encode(input);
    std::println("{}", encoded);

    auto decoded = decode(encoded);
    assert(input.size() == decoded.size());
    for (usize i = 0; i < input.size(); ++i) {
        assert(input[i] == decoded[i]);
    }
}
} // namespace Problem6

std::vector<int> productExceptSelf(std::vector<int> &nums) {
    std::vector<int> left;
    left.resize(nums.size());
    std::vector<int> right;
    right.resize(nums.size());
    std::vector<int> out;
    out.resize(nums.size());

    left[0] = nums[0];
    right[right.size() - 1] = nums[right.size() - 1];
    for (size_t i{1zu}; i < nums.size(); ++i) {
        left[i] = left[i - 1] * nums[i];
        const size_t r_idx = right.size() - 1 - i;
        right[r_idx] = right[r_idx + 1] * nums[r_idx];
    }
    out[0] = right[1];
    out[out.size() - 1] = left[out.size() - 2];
    for (size_t i{1zu}; i < out.size() - 1; ++i) {
        out[i] = left[i - 1] * right[i + 1];
    }
    return out;
}

bool is_valid_sodoku(const std::vector<std::vector<char>> board) {
    for (const auto &row : board) {
        int bits = 0;
        for (const char c : row) {
            if (c != '.') {
                const int flag = (1 << (c - '1'));
                if (bits & flag) {
                    return false;
                }
                bits |= flag;
            }
        }
    }
    for (usize col_idx{0zu}; col_idx < 9; ++col_idx) {
        int bits = 0;
        for (usize row_idx{0zu}; row_idx < 9; ++row_idx) {
            const char c = board[row_idx][col_idx];
            if (c != '.') {
                const int flag = (1 << (c - '1'));
                if (bits & flag) {
                    return false;
                }
                bits |= flag;
            }
        }
    }
    for (usize box_row{0zu}; box_row < 3; ++box_row) {
        for (usize box_col{0zu}; box_col < 3; ++box_col) {
            int bits = 0;
            for (usize row{0zu}; row < 3; ++row) {
                for (usize col{0zu}; col < 3; ++col) {
                    const usize row_g = 3 * box_row + row;
                    const usize col_g = 3 * box_col + col;
                    const char c = board[row_g][col_g];
                    if (c != '.') {
                        const int flag = (1 << (c - '1'));
                        if (bits & flag) {
                            return false;
                        }
                        bits |= flag;
                    }
                }
            }
        }
    }
    return true;
}

bool valid_parenthesis(std::string s) {
    std::array<char, 500> my_stack{};
    usize stack_ptr{0zu};
    std::array<std::pair<char, char>, 3> pairs{};
    pairs[0] = {'(', ')'};
    pairs[1] = {'[', ']'};
    pairs[2] = {'{', '}'};

    for (const char c : s) {
        if (stack_ptr >= my_stack.size()) {
            return false;
        }
        // Push
        if (c == '(' || c == '[' || c == '{') {
            my_stack[stack_ptr++] = c;
            continue;
        }
        // Pop
        if (stack_ptr == 0) {
            return false;
        }
        for (const auto &[left, right] : pairs) {
            if (c == right && my_stack[stack_ptr - 1] != left) {
                return false;
            }
        }
        --stack_ptr;
    }
    return (stack_ptr == 0);
}

namespace ProblemMinimumStack {
class MinStack {
    using Value = int;
    using MinValue = int;

public:
    MinStack() {
        min_stack_.reserve(128);
    }

    void push(int val) {
        if (min_stack_.empty()) {
            min_stack_.push_back({val, val});
        } else {
            const auto &[old_val, old_min] = min_stack_.back();
            if (old_min > val) {
                min_stack_.push_back({val, val});
            } else {
                min_stack_.push_back({val, old_min});
            }
        }
    }

    void pop() {
        min_stack_.pop_back();
    }

    int top() {
        return min_stack_.back().first;
    }

    int getMin() {
        return min_stack_.back().second;
    }

private:
    std::vector<std::pair<Value, MinValue>> min_stack_{};
};
} // namespace ProblemMinimumStack

int eval_prn(const std::vector<std::string> &tokens) {
    std::vector<int> values;
    values.reserve(tokens.size());

    for (usize i{0zu}; i < tokens.size(); ++i) {
        const std::string &t = tokens[i];
        if (t == "+" || t == "-" || t == "*" || t == "/") {
            const int right = values.back();
            values.pop_back();
            const int left = values.back();
            values.pop_back();
            if (t == "+") {
                values.push_back(left + right);
            } else if (t == "-") {
                values.push_back(left - right);
            } else if (t == "*") {
                values.push_back(left * right);
            } else if (t == "/") {
                values.push_back(left / right);
            }
        } else {
            values.push_back(std::stoi(t));
        }
    }
    return values.back();
}

std::vector<int> daily_temperature(std::vector<int> &temperatures) {
    using Position = uint32_t; // not using usize as 12 byte is awkward alignment
    std::vector<std::pair<int, Position>> my_stack;
    std::vector<int> out;
    out.resize(temperatures.size());

    my_stack.reserve(temperatures.size());
    for (usize i{0}; i < temperatures.size(); ++i) {
        const int value{temperatures[i]};
        while (!my_stack.empty() && value > my_stack.back().first) {
            const Position idx = my_stack.back().second;
            out[static_cast<usize>(idx)] = static_cast<int>(i - idx);
            my_stack.pop_back();
        }
        my_stack.push_back({value, static_cast<Position>(i)});
    }

    return out;
}

namespace ProblemIsPalindrome {
bool same_alphanum(char a, char b) {
    if (a >= 'A' && a <= 'Z')
        a += 'a' - 'A';
    if (b >= 'A' && b <= 'Z')
        b += 'a' - 'A';

    return a == b;
}

bool isPalindrome(const std::string &s) {
    usize left{0zu};
    usize right{s.size() - 1};

    while (left < right) {
        while (!std::isalnum(s[left])) {
            ++left;
        }
        while (!std::isalnum(s[right])) {
            if (right == 0) {
                return (left >= s.size());
            }
            --right;
        }
        if (left >= right) {
            return same_alphanum(s[left], s[right]);
        }
        if (!same_alphanum(s[left], s[right])) {
            return false;
        }
        ++left;
        if (right > 0) {
            --right;
        }
    }
    return true;
}

TreeNode *invert_tree(TreeNode *root) {
    if (!root)
        return nullptr;

    std::swap(root->left, root->right);

    invert_tree(root->left);
    invert_tree(root->right);

    return root;
}

} // namespace ProblemIsPalindrome

std::vector<int> two_sum2(const std::vector<int> &numbers, int target) {
    usize left{0zu};
    usize right{numbers.size() - 1};

    while (left < right) {
        while (left < right && numbers[left] + numbers[right] < target) {
            ++left;
        }
        if (numbers[left] + numbers[right] == target) {
            return {static_cast<int>(left + 1), static_cast<int>(right + 1)};
        }
        --right;
    }
    return {0, 0};
}

void rotate(std::vector<std::vector<int>> &matrix) {
    const usize n = matrix.size(); // assuming square matrix
    for (usize row{0zu}; row < n / 2; ++row) {
        for (usize col{row}; col < n - 1 - row; ++col) {
            const int tl = matrix[row][col];
            const int tr = matrix[col][n - 1 - row];
            const int br = matrix[n - 1 - row][n - 1 - col];
            const int bl = matrix[n - 1 - col][row];

            matrix[row][col] = bl;
            matrix[n - 1 - col][row] = br;
            matrix[n - 1 - row][n - 1 - col] = tr;
            matrix[col][n - 1 - row] = tl;
        }
    }
}

std::vector<int> plusOne(std::vector<int> &digits) {
    bool carry{true};              // The +1 is a type of carry
    std::vector<int> out = digits; // copy
    usize i = digits.size() - 1;
    while (true) {
        if (!carry) {
            break;
        }
        const int new_val = digits[i] + (carry ? 1 : 0);
        if (new_val == 10) {
            out[i] = 0;
        } else {
            out[i] = new_val;
            carry = false;
        }
        if (i == 0) {
            break;
        }
        --i;
    }
    if (carry) {
        out.insert(out.begin(), 1);
    }
    return out;
}

bool isSameTree(TreeNode *p, TreeNode *q) {
    if (p == nullptr) {
        return q == nullptr;
    }
    if (q == nullptr) {
        return p == nullptr;
    }
    return (p->val == q->val) && isSameTree(p->left, q->left) && isSameTree(p->right, q->right);
}

namespace ProblemBSTIsValid {
bool is_valid(TreeNode *node, int lower_bound, int upper_bound) {
    if (node == nullptr) {
        return true;
    }
    if (node->val <= lower_bound || node->val >= upper_bound) {
        return false;
    }
    return is_valid(
               node->left,
               lower_bound,
               node->val) &&
           is_valid(
               node->right,
               node->val,
               upper_bound);
}

bool isValidBST(TreeNode *root) {
    if (root == nullptr) {
        return true;
    }
    return is_valid(
               root->left,
               std::numeric_limits<int>::min(),
               root->val) &&
           is_valid(
               root->right,
               root->val,
               std::numeric_limits<int>::max());
}
} // namespace ProblemBSTIsValid

int singleNumber(std::vector<int> &nums) {
    int a = 0;
    for (const int i : nums) {
        a ^= i;
    }
    return a;
}

int getSum(int a, int b) {
    bool carry{false};
    int out{0};

    const unsigned int ua = static_cast<unsigned int>(a);
    const unsigned int ub = static_cast<unsigned int>(b);

    for (int i = 0; i < 32; ++i) {
        const unsigned int shifted = 1u << i;

        const bool abit = (ua & shifted) != 0;
        const bool bbit = (ub & shifted) != 0;

        if (abit || bbit) {     // at least one 1
            if (abit && bbit) { // both are 1
                if (carry) {    // (1 + 1), carry -> 1, carry
                    out |= shifted;
                }
                carry = true;
            } else {          // only one is 1
                if (!carry) { // (1 + 0), no carry -> 1, no carry
                    out |= shifted;
                    carry = false;
                } else { // (1 + 0), carry -> 0, carry
                    carry = true;
                }
            }
        } else if (carry) {
            out |= shifted;
            carry = false;
        }
    }
    return out;
}

struct TrieNode {
    std::vector<char> edges_{};
    std::vector<std::unique_ptr<TrieNode>> children_{};

    bool is_terminal_{false};

    void insert(std::string_view word) {
        if (word.empty()) {
            is_terminal_ = true;
            return;
        }

        for (usize i = 0; i < edges_.size(); ++i) {
            if (edges_[i] == word[0]) {
                assert(children_[i] && "Existing edges must contain children");
                children_[i]->insert(word.substr(1));
                return;
            }
        }
        edges_.push_back(word[0]);
        auto child = std::make_unique<TrieNode>();
        child->insert(word.substr(1));
        children_.push_back(std::move(child));
    }

    bool search(std::string_view word) const {
        if (word.empty()) {
            return is_terminal_;
        }
        for (usize i = 0; i < edges_.size(); ++i) {
            if (edges_[i] == word[0]) {
                assert(children_[i] && "Existing edges must contain children");
                return children_[i]->search(word.substr(1));
            }
        }
        return false;
    }

    bool starts_with(std::string_view word) const {
        if (word.empty()) {
            return true;
        }
        for (usize i = 0; i < edges_.size(); ++i) {
            if (edges_[i] == word[0]) {
                assert(children_[i] && "Existing edges must contain children");
                return children_[i]->starts_with(word.substr(1));
            }
        }
        return false;
    }

    void print_indent(int indent) const {
        for (int i = 0; i < indent; ++i) {
            std::print(" ");
        }
    }

    void print_edges(int indent = 0) const {
        for (usize i{0zu}; i < edges_.size(); ++i) {
            print_indent(indent);
            const char edge = edges_[i];
            std::println("{}", edge);
            const auto &child = children_[i];
            child->print_edges(indent + 4);
        }
    }

    bool is_terminal() const {
        return is_terminal_;
    }
};

void run_tests() {
    using namespace ds_neetcode;

    std::println("Running the tests.");
    if constexpr (test_problem_1) {
        problem_1_test();
    }
    if constexpr (test_problem_2) {
        problem_2_test();
    }
    if constexpr (test_problem_3) {
        problem_3_test();
    }
    if constexpr (test_problem_4) {
        problem_4_test();
    }
    if constexpr (test_problem_5) {
        problem_5_test();
    }
    if constexpr (test_problem_6) {
        using namespace Problem6;
        problem_6_test();
    }
    std::println("Passed all tests!");
}

std::vector<std::vector<int>> subsets(std::vector<int> &nums) {
    if (nums.empty()) {
        return {{}};
    }
    std::vector<int> new_vec{};
    new_vec.reserve(nums.size() - 1);
    for (usize i = 1; i < nums.size(); ++i) {
        new_vec.push_back(nums[i]);
    }

    auto ret = subsets(new_vec);

    const usize n = ret.size();
    for (usize i{0zu}; i < n; ++i) {
        std::vector<int> set2 = ret[i];
        set2.push_back(nums[0]);
        ret.push_back(set2);
    }

    return ret;
}

int maxProfit(std::vector<int> &prices) {
    int current_max{0};
    int current_accum{0};
    for (usize i{1zu}; i < prices.size(); ++i) {
        const int delta{prices[i] - prices[i - 1]};
        current_accum = std::max(0, current_accum + delta);
        current_max = std::max(current_max, current_accum);
    }
    return current_max;
}

namespace ProblemClimbingStairs {
#include <array>

class Solution {
public:
    int climbStairs(int n) {
        if (!inited_) {
            memo_.fill(-1);
            memo_[0] = 1;
            memo_[1] = 1;
            inited_ = true;
        }
        const usize idx = static_cast<usize>(n);
        if (memo_[idx] != -1) {
            return memo_[idx];
        }
        memo_[idx] = climbStairs(n - 1) + climbStairs(n - 2);
        return memo_[idx];
    }

private:
    std::array<int, 30> memo_;
    bool inited_{false};
};
} // namespace ProblemClimbingStairs

namespace ProblemNumIsland {
class Solution {
public:
    int numIslands(std::vector<std::vector<char>> &grid) {
        using usize = std::size_t;
        const usize n = grid.size();
        const usize m = grid[0].size();

        std::vector<std::vector<bool>> visited(n, std::vector<bool>(m, false));
        int count{0};

        for (usize row = 0; row < n; ++row) {
            for (usize col = 0; col < m; ++col) {
                if (grid[row][col] == '0') {
                    visited[row][col] = true;
                } else if (!visited[row][col]) { // Flood Fill
                    std::vector<std::pair<usize, usize>> to_check{};
                    to_check.push_back({row, col});
                    visited[row][col] = true;

                    while (!to_check.empty()) {
                        const auto [row_, col_] = to_check.back();
                        to_check.pop_back();

                        if (row_ > 0 && !visited[row_ - 1][col_] && grid[row_ - 1][col_] == '1') {
                            visited[row_ - 1][col_] = true;
                            to_check.push_back({row_ - 1, col_});
                        }
                        if (row_ < n - 1 && !visited[row_ + 1][col_] && grid[row_ + 1][col_] == '1') {
                            visited[row_ + 1][col_] = true;
                            to_check.push_back({row_ + 1, col_});
                        }
                        if (col_ > 0 && !visited[row_][col_ - 1] && grid[row_][col_ - 1] == '1') {
                            visited[row_][col_ - 1] = true;
                            to_check.push_back({row_, col_ - 1});
                        }
                        if (col_ < m - 1 && !visited[row_][col_ + 1] && grid[row_][col_ + 1] == '1') {
                            visited[row_][col_ + 1] = true;
                            to_check.push_back({row_, col_ + 1});
                        }
                    }

                    ++count;
                }
            }
        }

        return count;
    }
};
} // namespace ProblemNumIsland

using usize = std::size_t;
using isize = std::ptrdiff_t;

class DynamicArray {
public:
    DynamicArray(int capacity) {
        start_ = static_cast<int *>(std::malloc(static_cast<usize>(capacity) * sizeof(int)));
        end_ = start_;
        capacity_ = start_ + capacity;
    }

    ~DynamicArray() {
        std::free(start_);
    }

    int get(int i) const {
        return *(start_ + i);
    }

    void set(int i, int n) {
        *(start_ + i) = n;
    }

    void push_back(int n) {
        if (end_ == capacity_) {
            resize();
        }
        *end_ = n;
        ++end_;
    }

    int pop_back() {
        --end_;
        int val = *end_;
        return val;
    }

    void resize() {
        const int n_elements = get_size();
        const usize old_capacity = static_cast<usize>(capacity_ - start_);
        const usize new_capacity = static_cast<usize>(2 * static_cast<double>(old_capacity));

        auto tmp = static_cast<int *>(std::realloc(start_, new_capacity * sizeof(int)));
        assert(tmp);
        start_ = tmp;
        end_ = start_ + n_elements;
        capacity_ = start_ + new_capacity;
    }

    int get_size() const {
        return static_cast<int>(end_ - start_);
    }

    int get_capacity() const {
        return static_cast<int>(capacity_ - start_);
    }

    bool empty() const {
        return start_ == end_;
    }

private:
    int *start_{};
    int *end_{};
    int *capacity_{};
};

namespace ProblemLinkedListMerging {

struct ListNode {
    int val;
    ListNode *next;
    ListNode() : val(0), next(nullptr) {}
    ListNode(int x) : val(x), next(nullptr) {}
    ListNode(int x, ListNode *next) : val(x), next(next) {}
};

class Solution {
public:
    ListNode *mergeTwoLists(ListNode *list1, ListNode *list2) {
        if (!list1)
            return list2;
        if (!list2)
            return list1;

        ListNode *head = nullptr;
        ListNode *curr = nullptr;

        if (list1->val <= list2->val) {
            head = curr = list1;
            list1 = list1->next;
        } else {
            head = curr = list2;
            list2 = list2->next;
        }

        while (list1 && list2) {
            if (list1->val <= list2->val) {
                curr->next = list1;
                list1 = list1->next;
            } else {
                curr->next = list2;
                list2 = list2->next;
            }
            curr = curr->next;
        }
        while (list1) {
            curr->next = list1;
            list1 = list1->next;
            curr = curr->next;
        }
        while (list2) {
            curr->next = list2;
            list2 = list2->next;
            curr = curr->next;
        }
        return head;
    }
};
} // namespace ProblemLinkedListMerging

namespace ThreeSum {
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using T = int;
using usize = std::size_t;
using Triple = std::array<T, 3>;

constexpr bool k_print{false};

struct TripleHash {
    usize operator()(const Triple &t) const noexcept {
        auto hash = 0zu;
        hash = hash * 31 + t[0];
        hash = hash * 31 + t[1];
        hash = hash * 31 + t[2];
        return hash;
    }
};

class Solution {
public:
    vector<vector<T>> threeSum(vector<T> &nums) {
        static_assert(std::is_trivial_v<T>);
        auto n = nums.size();
        std::unordered_set<Triple, TripleHash> seen{};
        seen.reserve(4 * n);

        std::unordered_map<T, int> counts{};
        counts.reserve(n);
        for (T num : nums) {
            ++counts[num];
        }

        vector<vector<T>> out{};
        out.reserve(n);

        for (auto i{0zu}; i < n; ++i) {
            for (auto j{i + 1zu}; j < n; ++j) {
                auto x = nums[i];
                auto y = nums[j];
                auto target = -x - y;

                // need_x = 1 if x, y, target are distinct
                // need_x = 2 if ((x == y) || (x == target)) && (y == target)
                //            or equivalents (x == y) ^ (x == target)
                // need_y = 3 if x == y == target
                auto need_x = 1 + (x == y) + (x == target);

                // The y == x case is already covered before so we only need
                // to consider the case where (y != x) && (y == target);
                auto need_y = 1 + (y == target);

                // Checking for == 1 avoids having to do hashmap lookup
                auto enough_x = (need_x == 1) || (need_x <= counts[x]);
                auto enough_y = (need_y == 1) || (need_y <= counts[y]);
                if (!enough_x || !enough_y) {
                    continue;
                }

                auto it = counts.find(target);
                if (it == counts.end())
                    continue;
                assert(target == it->first);
                Triple t{x, y, target};

                std::sort(std::begin(t), std::end(t));
                if constexpr (k_print) {
                    std::cout << "t = (" << x << ", " << y << ", " << target << ")\n";
                }

                if (seen.find(t) == seen.end()) {
                    out.push_back({t[0], t[1], t[2]});
                    seen.insert(t);
                }
            }
        }
        return out;
    }
};
} // namespace ThreeSum
} // namespace ds_neetcode
