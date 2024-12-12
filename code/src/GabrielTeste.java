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
        double[][] matrix = readCSVToMatrix(csvPath);

        EigenDecomposition eigenDecomposition = decomposeMatrix(matrix);
        RealMatrix eigenVectors = eigenDecomposition.getV();
        RealMatrix eigenValues = eigenDecomposition.getD();
        double[][] eigenValuesArray = eigenValues.getData();
        int k = 1;
        getMaxAbsoluteValues(eigenValuesArray, k); //para buscar maior valor absoluto
        RealMatrix eigenVectorsTranspose = eigenDecomposition.getVT();
        RealMatrix originalMatrix = eigenVectors.multiply(eigenValues).multiply(eigenVectorsTranspose);

        RealMatrix matrix2 = new Array2DRowRealMatrix(matrix);

        printMatrix(matrix2, "Matrix (A)");
        printMatrix(eigenVectors, "EigenVectors (P)");
        printMatrix(eigenValues, "EigenValues (D)");
        printMatrix(eigenVectorsTranspose, "EigenVectorsTranspose (P^t)");
        printMatrix(originalMatrix, "OriginalMatrix (A = P.D.P^t)");


    }

    private static double getMaxAbsoluteValues(double[][] eigenValuesArray, int k) {

        double[] quantity = new double[k];

//        for (int i = 0; i < quantity.length; i++) {
//            quantity[i] = Double.MIN_VALUE;
//        }
//
//        for (int i = 0; i < eigenValuesArray.length; i++) {
//            if (Math.abs(eigenValuesArray[i][i]) > Math.abs(quantity[i])) {
//                quantity[i] = eigenValuesArray[i][i];
//            }
//        }
// return erra do, rodar matriz maior na diagonal, se for maior do que matrix K, insere, roda a matrix K por completo por pode ter sempre um maior valorAbsoluto.
        return manyValues;
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