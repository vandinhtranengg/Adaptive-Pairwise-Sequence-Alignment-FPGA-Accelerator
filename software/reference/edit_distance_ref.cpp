#include "edit_distance_ref.hpp"

#include <algorithm>
#include <numeric>
#include <vector>

using std::min;
using std::string;
using std::vector;
using std::iota;

uint32_t edit_distance_ref(
    const string& sequence_a,
    const string& sequence_b)
{
    size_t m = sequence_a.size();
    size_t n = sequence_b.size();

    vector<uint32_t> previous(n + 1);
    vector<uint32_t> current(n + 1);

    iota(previous.begin(), previous.end(), 0U);

    for (size_t i = 1; i <= m; i++) {
        current[0] = i;

        for (size_t j = 1; j <= n; j++) {
            uint32_t substitution_cost =
                (sequence_a[i - 1] == sequence_b[j - 1]) ? 0 : 1;

            uint32_t deletion =
                previous[j] + 1;

            uint32_t insertion =
                current[j - 1] + 1;

            uint32_t substitution =
                previous[j - 1] + substitution_cost;

            current[j] = min({
                deletion,
                insertion,
                substitution
            });
        }

        previous.swap(current);
    }

    return previous[n];
}