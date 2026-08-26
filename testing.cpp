#include <iostream>
#include <iomanip>
#include "sources.h"

int main() {
    double gauss = Sources().gaussianPulse(0.5, 0.5, 0.1);
    std::cout << std::fixed << std::setprecision(6) << gauss;
    return 0;
}
