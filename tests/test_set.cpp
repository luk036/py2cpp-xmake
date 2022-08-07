#include <doctest/doctest.h>  // for ResultBuilder, TestCase, CHECK, TEST_CASE

#include <py2cpp/set.hpp>  // for set

TEST_CASE("Test set") {
    const auto S = py::set<int>{1, 3, 4, 5, 1};
    auto count = 0;
    for (const auto& _ : S) {
        static_assert(sizeof(_) >= 0, "make compiler happy");
        ++count;
    }
    CHECK(count == 4);
}
