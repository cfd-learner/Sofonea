double*** createDistributionFunction(int xNum, int yNum);
double** createLayer(int xNum, int yNum);
void freeMemory(double*** arr, int xNum, int yNum);
void freeMemory(double** arr, int xNum, int yNum);
void Copy(double*** dest, double*** source, int xNum, int yNum);