#include <iostream>

int main() {

    int N1, N2; N1 = 8; N2 = 5;
    int n1, n2; n1 = 0; n2 = 0;
    int j = 0;

    for(n1=0; n1 < N1; n1++)
        for(n2=0; n2 < N2; n2++) {

            j = n1 * N2 + n2;

            std::cout << j << "(" << n1 << "," << n2 << ")" << std::endl;

        }

    std::cout << std::endl;

    // inverse
    int np = N1 * N2;
    for(j=0;j < np; j++) {

        n1 = j/N2; // j/N2 is integer division
        n2 = j-N2*(j/N2); // j/N2 is integer division

        std::cout << j << "(" << n1 << "," << n2 << ")" << std::endl;

    }
    
    return 0;

}

