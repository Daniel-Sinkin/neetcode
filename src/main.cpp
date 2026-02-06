// src/main.cpp

#include "pch.hpp" // IWYU pragma: keep
#include <cstddef>
#include <print>

namespace ds_neetcode {
struct ListNode {
    int val;
    ListNode *next;
    ListNode() : val(0), next(nullptr) {}
    ListNode(int x) : val(x), next(nullptr) {}
    ListNode(int x, ListNode *next) : val(x), next(next) {}
};

void print_ll(ListNode &root) {
    ListNode *curr{&root};
    while (curr) {
        int val = curr->val;
        curr = curr->next;
        if (curr) {
            std::print("{} -> ", val);
        } else {
            std::print("{}", val);
        }
    }
    std::println();
}
} // namespace ds_neetcode

int main() {
    using namespace ds_neetcode;

    ListNode x3{3};
    ListNode x2{2, &x3};
    ListNode x1{1, &x2};
    ListNode root{0, &x1};
    print_ll(root);

    ListNode *curr{&root};
    ListNode *prev{};
    while (curr) {
        ListNode *next{curr->next};
        curr->next = prev;
        prev = curr;
        curr = next;
    }
    print_ll(*prev);
}
