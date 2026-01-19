// src/main.cpp

#include "pch.hpp" // IWYU pragma: keep
#include <cmath>
#include <cstddef>

namespace ds_neetcode {
using usize = std::size_t;
using isize = std::ptrdiff_t;

namespace ProblemClimbingStairs {
#include <array>

class Solution {
public:
    int climbStairs(int n) {
        if(!inited_) {
            memo_.fill(-1);
            memo_[0] = 1;
            memo_[1] = 1;
            inited_ = true;
        }
        const usize idx = static_cast<usize>(n);
        if(memo_[idx] != -1) {
            return memo_[idx];
        }
        memo_[idx] = climbStairs(n - 1) + climbStairs(n - 2);
        return memo_[idx];
    }
private:
    std::array<int, 30> memo_;
    bool inited_{false};
};
} // namespaec ProblemClimbingStairs

int reverse_int(int x) {
    std::string res{std::to_string(x)};
    std::reverse(res.begin() + ((x < 0) ? 1 : 0), res.end());
    try {
        return std::stoi(res);
    } catch (std::out_of_range) {
        return 0;
    }
}

namespace ProblemNumIsland {
class Solution {
public:
    int numIslands(std::vector<std::vector<char>>& grid) {
        using usize = std::size_t;
        const usize n = grid.size();
        const usize m = grid[0].size();
    
        std::vector<std::vector<bool>> visited(n, std::vector<bool>(m, false));
        int count{0};

        for(usize row = 0; row < n; ++row) {
            for(usize col = 0; col < m; ++col) {
                if(grid[row][col] == '0') {
                    visited[row][col] = true;
                } else if (!visited[row][col]) { // Flood Fill
                    std::vector<std::pair<usize, usize>> to_check{};
                    to_check.push_back({row, col});
                    visited[row][col] = true;

                    while(!to_check.empty()) {
                        const auto [row_, col_] = to_check.back();
                        to_check.pop_back();

                        if(row_ > 0 && !visited[row_ - 1][col_] && grid[row_ - 1][col_] == '1') {
                            visited[row_ - 1][col_] = true;
                            to_check.push_back({row_ - 1, col_});
                        }
                        if(row_ < n - 1 && !visited[row_ + 1][col_] && grid[row_ + 1][col_] == '1') {
                            visited[row_ + 1][col_] = true;
                            to_check.push_back({row_ + 1, col_});
                        }
                        if(col_ > 0 && !visited[row_][col_ - 1] && grid[row_][col_ - 1] == '1') {
                            visited[row_][col_ - 1] = true;
                            to_check.push_back({row_, col_ - 1});
                        }
                        if(col_ < m - 1 && !visited[row_][col_ + 1] && grid[row_][col_ + 1] == '1') {
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

class DynamicArray {
public:
    DynamicArray(int capacity) {
        start_ = static_cast<int*>(std::malloc(static_cast<usize>(capacity) * sizeof(int)));
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
        if(end_ == capacity_) {
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

        auto tmp = static_cast<int*>(std::realloc(start_, new_capacity * sizeof(int)));
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
    int* start_{};
    int* end_{};
    int* capacity_{};
};

int trapping_rain_water(const std::vector<int>& height) {
    std::vector<int> prefix;
    prefix.resize(height.size(), 0);
    std::vector<int> suffix;
    suffix.resize(height.size(), 0);
    { // Filling prefix and suffix arrays
        int prefix_max{0};
        int suffix_max{0};
        for(usize i_prefix{0zu}; i_prefix < height.size(); ++i_prefix) {
            const usize i_suffix = height.size() - i_prefix - 1;

            prefix_max = std::max(height[i_prefix], prefix_max);
            suffix_max = std::max(height[i_suffix], suffix_max);

            prefix[i_prefix] = prefix_max;
            suffix[i_suffix] = suffix_max;
        }
    }
    int final_water{0};
    // Reusing prefix for final water level to avoid additional memory usage:w
    for(usize i{0zu}; i < prefix.size(); ++i) {
        final_water += std::min(prefix[i], suffix[i]) - height[i];
    }
    return final_water;
}

} // namespace ns_neetnode;

int main() {
    using namespace ds_neetcode;
}
