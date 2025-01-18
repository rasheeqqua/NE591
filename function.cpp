//
// Created by Hasibul H. Rasheeq on 1/17/25.
//

#include <cmath>
#include <vector>

std::vector<double> evaluate(std::vector<double> input) {
  std::vector<double> output;
  for (int i = 0; i < input.size(); i++) {
    double value = exp(input[i]);
    output.push_back(value);
  }

  return output;
}
