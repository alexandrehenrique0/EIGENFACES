import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

public class GabrielTeste {
    public static void main(String[] args) {
        String csvPath = "Input/TesteFuncao1/csv/exampleInputFunc2.csv";
        double[][] arrayOriginal = readCSVToMatrix(csvPath);
        int k = 3;

        EigenDecomposition eigenDecomposition = decomposeMatrix(arrayOriginal);
        // primeiro buscamos as matrizes com get, pois eigenDecomposition é uma matriz real para depois transformar em array
        RealMatrix eigenVectors = eigenDecomposition.getV();
        double[][] eigenVectorsArray = eigenVectors.getData();

        RealMatrix eigenValues = eigenDecomposition.getD();
        double[][] eigenValuesArray = eigenValues.getData();

        RealMatrix eigenVectorsTranspose = eigenDecomposition.getVT();
        double[][] eigenVectorsTransposeArray =  eigenVectorsTranspose.getData();

        double[][] valuesAndIndexArray = getMaxAbsoluteValues(eigenValuesArray, k); //para buscar maior valor absoluto

        //methods K
        double[][] eigenVectorsKArray = vectorsK(eigenVectorsArray, valuesAndIndexArray);

        double[][] absoluteValuesKArray = diagonalMatrix(valuesAndIndexArray);

        double[][] vectorsKTransposeArray = transposed_Matrix(eigenVectorsKArray);

        double[][] A = multiply_Matrices(multiply_Matrices(eigenVectorsArray, eigenValuesArray), eigenVectorsTransposeArray);

        double[][] Ak = multiply_Matrices(multiply_Matrices(eigenVectorsKArray, absoluteValuesKArray), vectorsKTransposeArray);

        double erroMedio = calculateEAM(arrayOriginal, Ak);

        print_Matrix(arrayOriginal, "Matrix original");
        print_Matrix(eigenVectorsArray, "Matrix P");
        print_Matrix(eigenValuesArray, "Matrix D");
        print_Matrix(eigenVectorsTransposeArray, "Matrix P^t");
        print_Matrix(A, "Matrix A = P.D.P^t");
        print_Matrix(arrayOriginal, "Matrix original");
        print_Matrix(eigenVectorsKArray, "Matrix Pk");
        print_Matrix(absoluteValuesKArray, "Matrix Dk");
        print_Matrix(vectorsKTransposeArray, "Matrix Pk^t");
        print_Matrix(Ak, "Matrix Ak = Pk.Dk.Pk^t");
        System.out.printf("Erro médio: %.3f\n",erroMedio);

    }

    private static void print_Line(int length, String pattern) {
        for (int i = 0; i < length; i++) {
            System.out.print(pattern);
        }
        System.out.println();
    }

    public static void print_Matrix(double[][] matrixToPrint, String matrixName) {
        System.out.println("Matriz: " + matrixName + " ↓");
        print_Line(matrixToPrint[0].length, "____________");

        for (double[] row : matrixToPrint) {
            System.out.print("|");
            for (int i = 0; i < row.length; i++) {
                System.out.printf("%8.3f\t", row[i]);
                if (i == row.length - 1) {
                    System.out.print("|");
                }
            }
            System.out.println();
        }
        print_Line(matrixToPrint[0].length, "============");
        System.out.println();
    }

    public static double calculateEAM(double[][] A, double[][] Ak) {
        int M = A.length;
        int N = A[0].length;
        double erroAbsMed = 0;
        // Percorre cada elemento da matriz
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                erroAbsMed += Math.abs(A[i][j] - Ak[i][j]);
            }
        }

        // Calcula o erro médio
        return erroAbsMed / (M * N);
    }

    public static double[][] multiply_Matrices(double[][] matrizLeft, double[][] matrizRight) {
        double[][] matrizResultante = new double[matrizLeft.length][matrizRight[0].length];
        for (int i = 0; i < matrizLeft.length; i++) {
            for (int j = 0; j < matrizRight[0].length; j++) {
                for (int k = 0; k < matrizRight.length; k++) {
                    matrizResultante[i][j] += matrizLeft[i][k] * matrizRight[k][j];
                }
            }
        }
        return matrizResultante;
    }

    public static double[][] transposed_Matrix(double[][] matriz) {
        double[][] matrizTransposta = new double[matriz[0].length][matriz.length];
        for (int i = 0; i < matriz.length; i++) {
            for (int j = 0; j < matriz[0].length; j++) {
                matrizTransposta[j][i] = matriz[i][j];
            }
        }
        return matrizTransposta;
    }


    //                  arrayValuesK (0 ,3)
    private static double[][] vectorsK(double[][] eigenVectorsArray, double[][] arrayValuesK) {
        double[][] newVectorsK = new double[eigenVectorsArray.length][arrayValuesK.length];

        for (int i = 0; i < eigenVectorsArray[0].length; i++) {
            for (int j = 0; j < arrayValuesK.length; j++) {
                newVectorsK[i][j] = eigenVectorsArray[i][(int) arrayValuesK[j][1]];
            }
        }


        return newVectorsK;
    }

    public static double[][] diagonalMatrix(double[][] matrixvaluesK) {
        double[][] matrixvaluesKPrint = new double[matrixvaluesK.length][matrixvaluesK.length];
        for (int i = 0; i < matrixvaluesK.length; i++) {
            matrixvaluesKPrint[i][i] = matrixvaluesK[i][0];
        }
        return matrixvaluesKPrint;
    }

    private static double[][] getMaxAbsoluteValues(double[][] eigenValuesArray, int k) {
        double[][] valuesAndIndexArray = new double[k][2];

        for (int i = 0; i < valuesAndIndexArray.length; i++) {
            valuesAndIndexArray[i][0] = Double.MIN_VALUE;
        }

        for (int i = 0; i < eigenValuesArray.length; i++) {
            double absValue = Math.abs(eigenValuesArray[i][i]);
            for (int j = 0; j < valuesAndIndexArray.length; j++) {
                if (absValue > Math.abs(valuesAndIndexArray[j][0])) {
                    for (int l = valuesAndIndexArray.length - 1; l > j; l--) {
                        valuesAndIndexArray[l][0] = valuesAndIndexArray[l - 1][0];
                        valuesAndIndexArray[l][1] = valuesAndIndexArray[l - 1][1];
                    }
                    valuesAndIndexArray[j][0] = eigenValuesArray[i][i];
                    valuesAndIndexArray[j][1] = i;
                    break;
                }
            }
        }

        return valuesAndIndexArray;
    }

    public static double[][] readCSVToMatrix(String filePath) {
        ArrayList<double[]> rows = new ArrayList<>();
        try (Scanner scanner = new Scanner(new File(filePath))) {
            while (scanner.hasNextLine()) {
                String line = scanner.nextLine().trim();
                if (!line.isEmpty()) {
                    String[] values = line.split(",");
                    double[] row = new double[values.length];
                    for (int i = 0; i < values.length; i++) {
                        row[i] = Double.parseDouble(values[i].trim());
                    }
                    rows.add(row);
                }
            }
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Erro: Arquivo CSV não encontrado no caminho: " + filePath, e);
        } catch (NumberFormatException e) {
            throw new RuntimeException("Erro: Valores inválidos encontrados no arquivo CSV.", e);
        }

        // Converte a lista de linhas em uma matriz
        double[][] matrix = new double[rows.size()][];
        for (int i = 0; i < rows.size(); i++) {
            matrix[i] = rows.get(i);
        }
        return matrix;
    }

    //* ---------------------------------------
    //* Decomposição de matriz simétrica
    public static EigenDecomposition decomposeMatrix(double[][] arrayParaDecompor) {
        Array2DRowRealMatrix matrix = new Array2DRowRealMatrix(arrayParaDecompor);
        EigenDecomposition eigenDecomposition = new EigenDecomposition(matrix);

//        RealMatrix eigenVectors = eigenDecomposition.getV();
//        RealMatrix eigenValues = eigenDecomposition.getD();
//        RealMatrix eigenVectorsTranspose = eigenDecomposition.getVT();

        return eigenDecomposition;
    }


    public static void printMatrix(RealMatrix matrix, String matrixName) {
        double[][] matrixToPrint = matrix.getData();
        System.out.println("Matriz: " + matrixName + " ↓");

        // Linha superior
        for (int i = 0; i < matrixToPrint[0].length; i++) {
            System.out.print("____________");
        }

        System.out.println(); // pula uma linha (estética)

        for (int j = 0; j < matrixToPrint.length; j++) {
            for (int i = 0; i < matrixToPrint[0].length; i++) {
                if (i == 0) {
                    System.out.printf("|%8.3f\t", matrixToPrint[j][i]);
                } else if (i == matrixToPrint[0].length - 1) {
                    System.out.printf("%8.3f\t|", matrixToPrint[j][i]);
                } else {
                    System.out.printf("%8.3f\t", matrixToPrint[j][i]);
                }
            }
            System.out.println();
        }
        // Linha final
        for (int i = 0; i < matrixToPrint[0].length; i++) {
            System.out.print("============");
        }
        System.out.println(); // pula duas linha (estética)
        System.out.println();
    }


    //! Error Messages
    public static void error_Location_Not_Found() {
        System.out.println("Localização inválida ou ficheiro não existe");
        System.exit(1);
    }

    public static void error_Invalid_Option() {
        System.out.println("Opção inválida");
        System.exit(1);
    }
    //! Error Messages
}