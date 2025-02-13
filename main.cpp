#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "DiffusionClass.cpp"

int main() {
    DiffusionSolver solver;
    std::string filepath;

    std::cout << "Choose input file type:\n"
              << "1. Standard input\n"
              << "2. Rotational symmetry test\n"
              << "3. Reflective symmetry test\n"
              << "Enter choice (1-3): ";

    int choice;
    std::cin >> choice;

    switch(choice) {
    case 1:
        filepath = "../input.txt";
        break;
    case 2:
        filepath = "../rotational_symmetry_test.txt";
        break;
    case 3:
        filepath = "../reflective_symmetry_test.txt";
        break;
    default:
        std::cerr << "Invalid choice\n";
        return 1;
    }

    try {
        if (!solver.readInput(filepath)) {
            return 1;
        }
        auto solution = solver.solve();
        solver.writeOutput(solution, "output.txt");
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}