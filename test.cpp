#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

// Program to generate input file for the S64 quadrature with 20,000 cells
int main() {
    std::ofstream input_file("input.txt");
    if (!input_file.is_open()) {
        std::cerr << "Error: Unable to create input file." << std::endl;
        return 1;
    }

    // Parameters
    int N = 32;                // Number of angles (S64 has 32 angles per hemisphere)
    int I = 20000;             // Number of spatial cells
    double sigma_T = 1.0;      // Total cross section
    double sigma_S = 0.5;      // Scattering cross section
    double q = 1.0;            // Source strength
    double L = 10.0;           // Slab width
    double eps_bar = 1.0e-5;   // Convergence tolerance
    int k_bar = 100;           // Maximum iterations

    // Write parameters to file
    input_file << N << " " << I << std::endl;
    input_file << sigma_T << " " << sigma_S << std::endl;
    input_file << q << " " << L << std::endl;
    input_file << eps_bar << " " << k_bar << std::endl;

    // S64 quadrature weights and points
    // First column is weight, second column is cosine of angle (mu)
    std::vector<std::pair<double, double>> quadrature = {
        {0.0243454785045699, 0.0243502926634244},
        {0.0242877337207517, 0.0729931217877990},
        {0.0241723811174015, 0.1214628192961200},
        {0.0239996942982292, 0.1696444204239920},
        {0.0237700828574152, 0.2174236437400070},
        {0.0234840914081050, 0.2646871622087670},
        {0.0231423982906572, 0.3113228719902110},
        {0.0227458139637091, 0.3572201583376680},
        {0.0222952790818783, 0.4022701579639910},
        {0.0217918622646617, 0.4463660172534640},
        {0.0212367575618268, 0.4894031457070530},
        {0.0206312816213117, 0.5312794640198940},
        {0.0199768705663602, 0.5718956462026340},
        {0.0192750765893078, 0.6111553551723930},
        {0.0185275642701200, 0.6489654712546570},
        {0.0177361066284412, 0.6852363130542330},
        {0.0169025809185708, 0.7198818501716100},
        {0.0160289641774258, 0.7528199072605310},
        {0.0151173285362012, 0.7839723589433410},
        {0.0141698363071298, 0.8132653151227970},
        {0.0131887348575274, 0.8406292962525800},
        {0.0121763512843555, 0.8659993981540920},
        {0.0111350869041916, 0.8893154459951140},
        {0.0100674115767651, 0.9105221370785020},
        {0.0089758578878487, 0.9295691721319390},
        {0.0078630152380124, 0.9464113748584020},
        {0.0067315239483593, 0.9610087996520530},
        {0.0055840697300656, 0.9733268277899110},
        {0.0044233799131820, 0.9833362538846260},
        {0.0032522289844892, 0.9910133714767440},
        {0.0020735166302813, 0.9963401167719550},
        {0.0008916403608482, 0.9993050417357720}
    };

    // Write quadrature data
    for (const auto& point : quadrature) {
        input_file << std::setprecision(16) << point.first << " " << point.second << std::endl;
    }

    input_file.close();
    std::cout << "Test input file for S64 with 20,000 cells generated successfully." << std::endl;

    return 0;
}