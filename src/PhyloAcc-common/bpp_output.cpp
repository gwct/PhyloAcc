#include "bpp_output.h"

#include <algorithm>

namespace phyloacc
{
double MedianInPlace(std::vector<double>& values, std::size_t begin, std::size_t end)
{
    std::sort(values.begin() + begin, values.begin() + end);
    return values[(begin + end) / 2];
}

double MeanRange(const std::vector<double>& values, std::size_t begin, std::size_t end)
{
    double total = 0;
    for (std::size_t i = begin; i < end; ++i)
        total += values[i];
    return total / (end - begin);
}

std::vector<std::vector<int> > CountZStates(int node_count,
                                            std::size_t begin,
                                            std::size_t end,
                                            const std::vector<std::vector<int> >& trace_Z,
                                            const std::vector<bool>& missing,
                                            int missing_count_value)
{
    std::vector<std::vector<int> > countZ(node_count, std::vector<int>(4, 0));
    for (int s = 0; s < node_count; ++s)
    {
        for (std::size_t i = begin; i < end; ++i)
            countZ[s][trace_Z[i][s] + 1]++;

        if (missing[s])
            countZ[s][0] = missing_count_value;
    }
    return countZ;
}

void WriteInitSummaryRow(std::ofstream& out_Z,
                         int cc,
                         double n_rate,
                         double c_rate,
                         double g_rate,
                         double l_rate,
                         double l2_rate,
                         const std::vector<std::vector<int> >& countZ,
                         int denom)
{
    out_Z << cc << "\t" << n_rate << "\t" << c_rate << "\t" << g_rate << "\t" << l_rate << "\t" << l2_rate;
    for (std::size_t s = 0; s < countZ.size(); ++s)
    {
        for (int k = 0; k < 4; ++k)
            out_Z << "\t" << static_cast<double>(countZ[s][k]) / denom;
    }
    out_Z << std::endl;
}
}
