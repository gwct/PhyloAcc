#ifndef PHYLOACC_COMMON_BPP_OUTPUT_H
#define PHYLOACC_COMMON_BPP_OUTPUT_H

#include <fstream>
#include <vector>

namespace phyloacc
{
double MedianInPlace(std::vector<double>& values, std::size_t begin, std::size_t end);
double MeanRange(const std::vector<double>& values, std::size_t begin, std::size_t end);

std::vector<std::vector<int> > CountZStates(int node_count,
                                            std::size_t begin,
                                            std::size_t end,
                                            const std::vector<std::vector<int> >& trace_Z,
                                            const std::vector<bool>& missing,
                                            int missing_count_value);

void WriteInitSummaryRow(std::ofstream& out_Z,
                         int cc,
                         double n_rate,
                         double c_rate,
                         double g_rate,
                         double l_rate,
                         double l2_rate,
                         const std::vector<std::vector<int> >& countZ,
                         int denom);
}

#endif
