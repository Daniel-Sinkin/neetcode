// src/main.cpp
#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cstddef>
#include <print>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ds_neetcode {

using usize = std::size_t;

using lc_chars = std::array<int, 26>; // lower case ascii chars

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

Encoded encode(const Decoded& strs) {
    if(strs.empty()) {
        return "";
    }
    std::string out;
    for(const std::string& word : strs) {
        const std::string word_size_s = std::to_string(word.size());
        if(word_size_s.size() == 1) {
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

Decoded decode(const Encoded& s) {
    Decoded out;
    usize i{0zu};
    while(i < s.size()) {
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
    for(usize i = 0; i < input.size(); ++i) {
        assert(input[i] == decoded[i]);
    }
}
} // namespace Problem6
} // namespace ds_neetcode


int main() {
    using namespace ds_neetcode;

    if constexpr (true) {
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
}
