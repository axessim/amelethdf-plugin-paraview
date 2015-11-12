/*
 * Copyright 2012 <nathanael.muot@axessim.fr>
 *
 */


// Declares a test suite named "sanity_check"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(sanity_check)

// Declares a test case named "sanity_check"
BOOST_AUTO_TEST_CASE(sanity_check) {
    BOOST_CHECK(true);  // We certainly hope that true is true
    BOOST_CHECK_EQUAL(2, 1+1);  // The value 1+1 should equal 2

    int x[] = {1, 2, 3};
    int y[] = {1, 2, 3};
    /* These arrays of length 3 are equal */
    BOOST_CHECK_EQUAL_COLLECTIONS(x, x+3, y, y+3);

    double a = 1.51;
    double b = 1.52;
    BOOST_CHECK_CLOSE_FRACTION(a, b, 0.11);  // These equal within 0.1
}

BOOST_AUTO_TEST_SUITE_END()
