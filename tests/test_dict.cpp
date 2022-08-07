#include <doctest/doctest.h>  // for ResultBuilder, TestCase, CHECK, TEST_CASE

#include <py2cpp/dict.hpp>  // for dict, key_iterator
#include <unordered_map>    // for operator!=, __hash_map_const_iterator
#include <utility>          // for pair

TEST_CASE("Test set") {
    using E = std::pair<double, int>;
    const auto S = py::dict<double, int>{E{0.1, 1}, E{0.3, 3}, E{0.4, 4}};
    auto count = 0;
    for (const auto& p : S) {
        static_assert(sizeof(p) >= 0, "make compiler happy");
        ++count;
    }
    CHECK(count == 3);
}
