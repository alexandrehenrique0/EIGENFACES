import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import java.io.IOException;
import java.io.File;
import java.util.Arrays;

public class Andre_Alexandre_RitaTest {
    public static double[][] matrizBase7x7 = {
            {1, 2, 3, 4, 5, 6, 7},
            {8, 9, 10, 11, 12, 13, 14},
            {15, 16, 17, 18, 19, 20, 21},
            {22, 23, 24, 25, 26, 27, 28},
            {29, 30, 31, 32, 33, 34, 35},
            {36, 37, 38, 39, 40, 41, 42},
            {43, 44, 45, 46, 47, 48, 49}
    };
    public static void main(String[] args) {

        testarConvertImagesToVectors();
        testarBuildImageMatrixM();
        testarCalcularVetorMedio();
        testarCentralizarMatriz();
        testarCalcularCovariancia();
        testarNormalizeEigenVectors();
        testarPrintMatrix();
    }

    // Teste do método convertImagesToVectors
    public static void testarConvertImagesToVectors() {
        System.out.println("A testar convertImagesToVectors com 40 ficheiros...");

        // Caminho para a pasta contendo os 40 ficheiros CSV
        String csvPath = "C:/Users/andre/IdeaProjects/lapr1-24-25_DAB_02/Input/Funcao2-3/csv";

        try {
            // Chama o método para converter os ficheiros CSV em vetores
            int[][] result = Andre_Alexandre_Rita.convertImagesToVectors(csvPath);

            // Verificação 1: Se a matriz tem exatamente 40 linhas (40 ficheiros)
            if (result.length == 40) {
                System.out.println("convertImagesToVectors: Teste de quantidade de ficheiros bem sucedido!");

                // Verificação 2: Se cada linha (vetor) tem o mesmo número de colunas
                int numColunas = result[0].length; // Pega o tamanho do primeiro vetor
                boolean colunasConsistentes = true;

                for (int i = 0; i < result.length; i++) {
                    if (result[i].length != numColunas) {
                        colunasConsistentes = false;
                        System.out.println("Falha: Vetor no ficheiro " + (i + 1) + " tem tamanho inconsistente.");
                        break;
                    }
                }

                if (colunasConsistentes) {
                    System.out.println("convertImagesToVectors: Teste de consistência de colunas bem sucedido!");
                } else {
                    System.out.println("convertImagesToVectors: Falha - Inconsistência no tamanho dos vetores.");
                }

            } else {
                System.out.println("convertImagesToVectors: Falha - Esperava 40 ficheiros, mas obteve " + result.length + ".");
            }

        } catch (IOException e) {
            System.out.println("convertImagesToVectors: Falha - Exceção lançada: " + e.getMessage());
        }
    }


    // Teste do método buildImageMatrixM
    public static void testarBuildImageMatrixM() {
        System.out.println("A testar buildImageMatrixM...");
        int[][] imageVectors = {
                {1, 2, 3, 4, 5, 6, 7},
                {8, 9, 10, 11, 12, 13, 14},
                {15, 16, 17, 18, 19, 20, 21},
                {22, 23, 24, 25, 26, 27, 28},
                {29, 30, 31, 32, 33, 34, 35},
                {36, 37, 38, 39, 40, 41, 42},
                {43, 44, 45, 46, 47, 48, 49}
        };
        double[][] expected = {
                {1, 8, 15, 22, 29, 36, 43},
                {2, 9, 16, 23, 30, 37, 44},
                {3, 10, 17, 24, 31, 38, 45},
                {4, 11, 18, 25, 32, 39, 46},
                {5, 12, 19, 26, 33, 40, 47},
                {6, 13, 20, 27, 34, 41, 48},
                {7, 14, 21, 28, 35, 42, 49}
        };

        double[][] result = Andre_Alexandre_Rita.buildImageMatrixM(imageVectors);

        if (matrizesIguais(result, expected)) {
            System.out.println("buildImageMatrixM: Teste bem sucedido!");
        } else {
            System.out.println("buildImageMatrixM: Falha - Resultado incorreto.");
        }
    }

    // Teste do método calcularVetorMedio
    public static void testarCalcularVetorMedio() {
        System.out.println("A testar calcularVetorMedio...");
        double[] expected= {4.0, 11.0, 18.0, 25.0, 32.0, 39.0, 46.0}; // Média das linhas

        double[][] matrizBase7x7 = {
                {1, 2, 3, 4, 5, 6, 7},
                {8, 9, 10, 11, 12, 13, 14},
                {15, 16, 17, 18, 19, 20, 21},
                {22, 23, 24, 25, 26, 27, 28},
                {29, 30, 31, 32, 33, 34, 35},
                {36, 37, 38, 39, 40, 41, 42},
                {43, 44, 45, 46, 47, 48, 49}
        };
        double[] result = Andre_Alexandre_Rita.calcularVetorMedio(matrizBase7x7);

        if (vetoresIguais(result, expected)) {
            System.out.println("calcularVetorMedio: Teste bem sucedido!");
        } else {
            System.out.println("calcularVetorMedio: Falha - Resultado incorreto.");
        }
    }

    // Teste do método centralizarMatriz
    public static void testarCentralizarMatriz() {
        System.out.println("A testar centralizarMatriz...");
        double[] vetorMedio = {4.0, 11.0, 18.0, 25.0, 32.0, 39.0, 46.0};


        double[][] expected = {
                {-3, -2, -1, 0, 1, 2, 3},
                {-3, -2, -1, 0, 1, 2, 3},
                {-3, -2, -1, 0, 1, 2, 3},
                {-3, -2, -1, 0, 1, 2, 3},
                {-3, -2, -1, 0, 1, 2, 3},
                {-3, -2, -1, 0, 1, 2, 3},
                {-3, -2, -1, 0, 1, 2, 3}
        };

        double[][] result = Andre_Alexandre_Rita.centralizarMatriz(matrizBase7x7, vetorMedio);

        if (matrizesIguais(result, expected)) {
            System.out.println("centralizarMatriz: Teste bem sucedido!");
        } else {
            System.out.println("centralizarMatriz: Falha - Resultado incorreto.");
        }
    }

    // Teste do método calcularCovariancia
    public static void testarCalcularCovariancia() {
        System.out.println("A testar calcularCovariancia...");
        double[][] matrizCentralizada = {
                {-3, -2, -1, 0, 1, 2, 3},
                {-3, -2, -1, 0, 1, 2, 3},
                {-3, -2, -1, 0, 1, 2, 3},
                {-3, -2, -1, 0, 1, 2, 3},
                {-3, -2, -1, 0, 1, 2, 3},
                {-3, -2, -1, 0, 1, 2, 3},
                {-3, -2, -1, 0, 1, 2, 3}
        };
        double[][] expected = {
                {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0},
                {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0},
                {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0},
                {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0},
                {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0},
                {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0},
                {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0}
        };

        double[][] result = Andre_Alexandre_Rita.calcularCovariancia(matrizCentralizada);

        if (matrizesIguais(result, expected)) {
            System.out.println("calcularCovariancia: Teste bem sucedido!");
        } else {
            System.out.println("calcularCovariancia: Falha - Resultado incorreto.");
        }
    }

    // Teste do método normalizeEigenVectors
    public static void testarNormalizeEigenVectors() {
        System.out.println("A testar normalizeEigenVectors...");
        double[][] eigenVectors = {
                {3, 4},
                {0, 0}
        };

        double[][] expected = {
                {1.0, 1.0},
                {0.0, 0.0}
        };

        RealMatrix matrix = new org.apache.commons.math3.linear.Array2DRowRealMatrix(eigenVectors);
        Andre_Alexandre_Rita.normalizeEigenVectors(matrix);
        double[][] result = matrix.getData();

        if (matrizesIguais(result, expected)) {
            System.out.println("normalizeEigenVectors: Teste bem sucedido!");
        } else {
            System.out.println("normalizeEigenVectors: Falha - Resultado incorreto.");
        }
    }

    // Teste do método printMatrix
    public static void testarPrintMatrix() {
        System.out.println("A testar printMatrix...");
        double[][] matriz = {
                {1.0, 2.0},
                {3.0, 4.0}
        };

        System.out.println("Saída esperada:");
        System.out.println("[1.0, 2.0]");
        System.out.println("[3.0, 4.0]");
        System.out.println("Saída do método:");
        Andre_Alexandre_Rita.printMatrix(matriz, "Matriz imprimida:");
    }

    // Funções auxiliares para comparar matrizes e vetores
    private static boolean matrizesIguais(double[][] matriz1, double[][] matriz2) {
        if (matriz1.length != matriz2.length || matriz1[0].length != matriz2[0].length) {
            return false;
        }

        for (int i = 0; i < matriz1.length; i++) {
            for (int j = 0; j < matriz1[0].length; j++) {
                if (Math.abs(matriz1[i][j] - matriz2[i][j]) > 1e-6) {
                    return false;
                }
            }
        }
        return true;
    }

    private static boolean vetoresIguais(double[] vetor1, double[] vetor2) {
        if (vetor1.length != vetor2.length) {
            return false;
        }

        for (int i = 0; i < vetor1.length; i++) {
            if (Math.abs(vetor1[i] - vetor2[i]) > 1e-6) {
                return false;
            }
        }
        return true;
    }
}




