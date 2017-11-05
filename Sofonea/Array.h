
double*** createArray(int basisNum, int xNum, int yNum);
double** createArray(int xNum, int yNum);
void freeMemory(double*** arr, int basisNum, int xNum);
void freeMemory(double** arr, int xNum);
void freeMemory(double* arr);
void freeMemory(int* arr);
void Copy(double*** dest, double*** source, int basisNum, int xNum, int yNum);