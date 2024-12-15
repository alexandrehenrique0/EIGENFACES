import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

//! FUNCIONALIDADE 1, DECOMPOR UMA MATRIZ SIMETRICA
public class GabrielFuncionalidadeUm {
    public static void main(String[] args) {
        String csvPath = "Input/Funcao1/exampleInputFunc2.csv";
        double[][] arrayOriginal = readCSVToArray(csvPath);
        int k = 1;

        EigenDecomposition eigenDecomposition = decomposeMatrix(arrayOriginal);
        double[][] eigenVectorsArray = getEigenVectors(eigenDecomposition);
        double[][] eigenValuesArray = getEigenValues(eigenDecomposition);
        double[][] eigenVectorsTransposeArray = getEigenVectorsTranspose(eigenDecomposition);

        double[][] valuesAndIndexArray = getValuesAndIndexArray(eigenValuesArray, k); //para buscar maior valor absoluto

        //methods K
        double[][] eigenVectorsKArray = getEigenVectorsKArray(eigenVectorsArray, valuesAndIndexArray);

        double[][] absoluteValuesKArray = constructDiagonalMatrix(valuesAndIndexArray);

        double[][] vectorsKTransposeArray = transposed_Matrix(eigenVectorsKArray);

        double[][] A = multiplyVectorsValuesVectorsTransposed(eigenVectorsArray, eigenValuesArray, eigenVectorsTransposeArray);

        double[][] Ak = multiplyVectorsValuesVectorsTransposed(eigenVectorsKArray, absoluteValuesKArray, vectorsKTransposeArray);

        double erroMedio = calculateEAM(arrayOriginal, Ak);

        saveMatrixToFile(Ak, csvPath, "Output/Func1");


        print_Matrix(arrayOriginal, "Matrix original");
        print_Matrix(eigenVectorsArray, "Matrix P");
        print_Matrix(eigenValuesArray, "Matrix D");
        print_Matrix(eigenVectorsTransposeArray, "Matrix P^t");
        print_Matrix(A, "Matrix A = P.D.P^t");
        print_Matrix(arrayOriginal, "Matrix original");
        print_Matrix(eigenVectorsKArray, "Matrix Pk");
        print_Matrix(absoluteValuesKArray, "Matrix Dk");
        print_Matrix(valuesAndIndexArray, "Valores absolutos e índices");
        print_Matrix(vectorsKTransposeArray, "Matrix Pk^t");
        print_Matrix(Ak, "Matrix Ak = Pk.Dk.Pk^t");
        System.out.printf("Erro médio: %.3f\n", erroMedio);

    }

    private static void saveMatrixToFile(double[][] matrix, String inputCsvPath, String outputFolderPath) {
        File outputFolder = new File(outputFolderPath);
        if (!outputFolder.exists()) {
            if (outputFolder.mkdirs()) {
                System.out.println("Diretório criado: " + outputFolderPath);
            } else {
                System.err.println("Falha ao criar o diretório: " + outputFolderPath);
                return;
            }
        }

        String csvFileName = new File(inputCsvPath).getName();
        String newFileName = "Reconstruida-" + csvFileName;
        String outputPath = outputFolderPath + "/" + newFileName;

        int counter = 1;
        File file = new File(outputPath);
        while (file.exists()) {
            file = new File(outputFolderPath + "/" + newFileName.replace(".csv", "(" + counter + ").csv"));
            counter++;
        }

        try (PrintWriter writer = new PrintWriter(file)) {
            for (double[] row : matrix) {
                String rowString = String.join(" ; ", Arrays.stream(row)
                        .mapToObj(val -> String.format("%.3f", val))
                        .toArray(String[]::new));
                writer.println(rowString);
            }
            System.out.println("Arquivo CSV criado com sucesso: " + file.getName());
        } catch (IOException e) {
            System.err.println("Erro ao salvar a matriz no arquivo: " + e.getMessage());
            e.printStackTrace();
        }
    }

    // serve para A e para Ak
    public static double[][] multiplyVectorsValuesVectorsTransposed(double[][] matrixVectors, double[][] matrixValues, double[][] matrixVectorsTranspose) {
        return multiply_Matrices(multiply_Matrices(matrixVectors, matrixValues), matrixVectorsTranspose);
    }

    public static double[][] getEigenVectors(EigenDecomposition eigenDecomposition) {
        RealMatrix eigenVectors = eigenDecomposition.getV();
        return eigenVectors.getData();
    }

    public static double[][] getEigenValues(EigenDecomposition eigenDecomposition) {
        RealMatrix eigenValues = eigenDecomposition.getD();
        return eigenValues.getData();
    }

    public static double[][] getEigenVectorsTranspose(EigenDecomposition eigenDecomposition) {
        RealMatrix eigenVectorsTranspose = eigenDecomposition.getVT();
        return eigenVectorsTranspose.getData();
    }

    private static void print_Line(int length, String pattern) {
        for (int i = 0; i < length; i++) {
            System.out.print(pattern);
        }
        System.out.println();
    }

    private static void print_Matrix(double[][] matrixToPrint, String matrixName) {
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

    private static double calculateEAM(double[][] A, double[][] Ak) {
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

    private static double[][] multiply_Matrices(double[][] matrizLeft, double[][] matrizRight) {
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

    private static double[][] transposed_Matrix(double[][] matriz) {
        double[][] matrizTransposta = new double[matriz[0].length][matriz.length];
        for (int i = 0; i < matriz.length; i++) {
            for (int j = 0; j < matriz[0].length; j++) {
                matrizTransposta[j][i] = matriz[i][j];
            }
        }
        return matrizTransposta;
    }

    private static double[][] getEigenVectorsKArray(double[][] eigenVectorsArray, double[][] arrayValuesK) {
        double[][] newVectorsK = new double[eigenVectorsArray.length][arrayValuesK.length];
        //TODO ao inves de dar system exit, pode retornar a pergunta de quantos K o utilizador quer
        if (arrayValuesK.length > eigenVectorsArray[0].length) {
            System.out.println("O valor de k não pode ser maior que o número de colunas da matriz de autovetores.");
            System.exit(1);
        } else if (arrayValuesK.length < 1) {
            System.out.println("O valor de k não pode ser menor que 1.");
            System.exit(1);
        }

        for (int i = 0; i < eigenVectorsArray[0].length; i++) {
            for (int j = 0; j < arrayValuesK.length; j++) {
                newVectorsK[i][j] = eigenVectorsArray[i][(int) arrayValuesK[j][1]];
            }
        }


        return newVectorsK;
    }

    private static double[][] constructDiagonalMatrix(double[][] matrixvaluesK) {
        double[][] matrixvaluesKPrint = new double[matrixvaluesK.length][matrixvaluesK.length];
        for (int i = 0; i < matrixvaluesK.length; i++) {
            matrixvaluesKPrint[i][i] = matrixvaluesK[i][0];
        }
        return matrixvaluesKPrint;
    }

    private static double[][] getValuesAndIndexArray(double[][] eigenValuesArray, int k) {
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

    private static double[][] readCSVToArray(String filePath) {
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
    private static EigenDecomposition decomposeMatrix(double[][] arrayParaDecompor) {
        Array2DRowRealMatrix matrix = new Array2DRowRealMatrix(arrayParaDecompor);
        EigenDecomposition eigenDecomposition = new EigenDecomposition(matrix);

        return eigenDecomposition;
    }
}


//        System.out.println(averageColumn.length);
//        print_Matrix(deviationMatrixA, "Matriz de Desvios (matrixTotalImages)");
//        print_Matrix(deviationMatrixATranspose, "Matriz Transposta de Desvios (matrixTotalImages^T)");
//        print_Matrix(covariance, "Matriz de Covariância (C)");
//        print_Matrix(eigenValuesATxA, "Valores Próprios de matrixTotalImages^T . matrixTotalImages");
//        print_Matrix(eigenVectorsATxA, "Vetores Próprios de matrixTotalImages^T . matrixTotalImages");
//        print_Matrix(eigenVectorsAxAT, "Vetores Próprios de matrixTotalImages . matrixTotalImages^T");
//        print_Matrix(eigenValuesAxAT, "Valores Próprios de matrixTotalImages . matrixTotalImages^T");
//        print_Matrix(normalizedEigenVectors, "Vetores Próprios Normalizados");
//        System.out.println("Tamanho do linearizedOriginalMatrix: " + linearizedOriginalMatrix.length);
//        System.out.println("Tamanho do vetor médio (meanVector): " + averageColumn.length);
//        System.out.println("Tamanho de phi ϕ : " + phi.length);
//        System.out.println("Pesos Calculados (W): " + Arrays.toString(weights));
//        System.out.println("Tamanho de reconstructedImageArray: " + reconstructedImageArray.length);
//        print_Matrix(new double[][]{reconstructedImageArray}, "Imagem Reconstruída (1D)");
//        System.out.println("Tamanho de reconstructedImageArray: " + reconstructedImageArray.length);
//        System.out.println("Tamanho esperado: " + (originalArray.length * originalArray[0].length));
//        System.out.println("Dimensões da imagem original: " + originalArray.length + "x" + originalArray[0].length);
//        System.out.println("Tamanho esperado: " + (originalArray.length * originalArray[0].length));
//        System.out.println("Tamanho do vetor reconstruído: " + reconstructedImageArray.length);
//        System.out.println("Tamanho do vetor médio (meanVector): " + averageColumn.length);
//        System.out.println("Dimensões de normalizedVectors: " + normalizedEigenVectors.length + "x" + normalizedEigenVectors[0].length);
//        System.out.println("Tamanho dos pesos (weights): " + weights.length);
//        print_Matrix(reconstructedImageMatrix, "Imagem Reconstruída (2D)");