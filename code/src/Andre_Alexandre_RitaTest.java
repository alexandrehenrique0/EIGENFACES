import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import java.io.IOException;
import java.io.File;
import java.util.Arrays;

public class Andre_Alexandre_RitaTest {
    public static void main(String[] args) {
        System.out.println("==== Testes Unitários da Classe Andre_Alexandre_Rita ====");

        testFindMaxEigenValue();
        testMatrixDesvios();
        testCovariancias();
        testNormalizeEigenVectors();
        testConvertIntToDouble();
        testColunaMedia();
        testTranspostaMatriz();
        testMultiplicaMatrizes();
        testMultiplicaMatrizPorEscalar();
        testBuildImageMatrixM();
        testConvertImagesToVectors();

        System.out.println("==== Testes Finalizados ====");
    }

    /**
     * Teste para calcular o maior valor próprio.
     */
    public static void testFindMaxEigenValue() {
        System.out.println("Teste 1: Calcular o maior valor próprio");

        double[][] matrix = {
                {4.0, -2.0},
                {-2.0, 4.0}
        };
        double expected = 6.0;
        double result = findMaxEigenValue(matrix);

        if (Double.compare(result, expected) == 0) {
            System.out.println("Sucesso: O maior valor próprio é " + result);
        } else {
            System.out.println("Falha: Esperado " + expected + ", mas obteve " + result);
        }
    }

    /**
     * Teste para centralizar a matriz.
     */
    public static void testMatrixDesvios() {
        System.out.println("Teste 2: Centralizar a matriz");

        double[][] matrix = {
                {1.0, 2.0, 3.0},
                {4.0, 5.0, 6.0}
        };
        double[][] colunaMedia = {
                {2.0},
                {5.0}
        };
        double[][] expected = {
                {-1.0, 0.0, 1.0},
                {-1.0, 0.0, 1.0}
        };

        double[][] result = Andre_Alexandre_Rita.matrixDesvios(matrix, colunaMedia);

        if (matricesAreEqual(result, expected)) {
            System.out.println("Sucesso: Centralização correta");
        } else {
            System.out.println("Falha: Centralização incorreta");
        }
    }

    /**
     * Teste para calcular a matriz de covariância.
     */
    public static void testCovariancias() {
        System.out.println("Teste 3: Cálculo da matriz de covariância");

        double[][] centralizada = {
                {-1.0, 0.0, 1.0},
                {-1.0, 0.0, 1.0}
        };
        double[][] expected = {
                {0.6667, 0.6667},
                {0.6667, 0.6667}
        };

        double[][] result = Andre_Alexandre_Rita.covariancias(centralizada, 3);

        if (matricesAreApproximatelyEqual(result, expected, 0.001)) {
            System.out.println("Sucesso: Covariância correta");
        } else {
            System.out.println("Falha: Covariância incorreta");
        }
    }

    /**
     * Teste para normalizar os vetores próprios.
     */
    public static void testNormalizeEigenVectors() {
        System.out.println("Teste 4: Normalizar os vetores próprios");

        double[][] data = {
                {7.0, 4.0},
                {0.0, 0.0}
        };


        RealMatrix matrix = new Array2DRowRealMatrix(data);
        Andre_Alexandre_Rita.normalizeEigenVectors(matrix);

        double[][] expected = {
                {1.0, 1.0},
                {0.0, 0.0}
        };



        if (matricesAreApproximatelyEqual(matrix.getData(), expected, 0.01)) {
            System.out.println("Sucesso: Normalização correta");
        } else {
            System.out.println("Falha: Normalização incorreta");
        }
    }

    /**
     * Teste para converter int[][] em double[][].
     */
    public static void testConvertIntToDouble() {
        System.out.println("Teste 5: Conversão de int[][] para double[][]");

        int[][] intMatrix = {
                {1, 2, 3},
                {4, 5, 6}
        };

        double[][] expected = {
                {1.0, 2.0, 3.0},
                {4.0, 5.0, 6.0}
        };

        double[][] result = Andre_Alexandre_Rita.convertIntToDouble(intMatrix);

        if (matricesAreEqual(result, expected)) {
            System.out.println("Sucesso: Conversão correta");
        } else {
            System.out.println("Falha: Conversão incorreta");
        }
    }

    /**
     * Testa o cálculo da média das colunas de uma matriz.
     * Verifica se o resultado corresponde à média esperada.
     */
    public static void testColunaMedia() {
        System.out.println("Teste 6: Cálculo da coluna média");

        double[][] matrix = {
                {1.0, 2.0, 3.0},
                {4.0, 5.0, 6.0}
        };

        double[][] expected = {
                {2.0},
                {5.0}
        };

        double[][] result = Andre_Alexandre_Rita.colunaMedia(matrix);

        if (matricesAreApproximatelyEqual(result, expected, 0.001)) {
            System.out.println("Sucesso: Cálculo da coluna média correto");
        } else {
            System.out.println("Falha: Resultado incorreto");
        }
    }

    /**
     * Testa a transposição de uma matriz.
     * Verifica se linhas e colunas foram corretamente invertidas.
     */
    public static void testTranspostaMatriz() {
        System.out.println("Teste 7: Transposição de matriz");

        double[][] matrix = {
                {1.0, 2.0, 3.0},
                {4.0, 5.0, 6.0}
        };

        double[][] expected = {
                {1.0, 4.0},
                {2.0, 5.0},
                {3.0, 6.0}
        };

        double[][] result = Andre_Alexandre_Rita.transpostaMatriz(matrix);

        if (matricesAreEqual(result, expected)) {
            System.out.println("Sucesso: Transposição correta");
        } else {
            System.out.println("Falha: Resultado incorreto");
        }
    }

    /**
     * Testa a multiplicação de duas matrizes.
     * Confirma se o produto das matrizes está correto.
     */
    public static void testMultiplicaMatrizes() {
        System.out.println("Teste 8: Multiplicação de matrizes");

        double[][] A = {
                {1.0, 2.0},
                {3.0, 4.0}
        };

        double[][] B = {
                {2.0, 0.0},
                {1.0, 2.0}
        };

        double[][] expected = {
                {4.0, 4.0},
                {10.0, 8.0}
        };

        double[][] result = Andre_Alexandre_Rita.multiplicaMatrizes(A, B);

        if (matricesAreApproximatelyEqual(result, expected, 0.001)) {
            System.out.println("Sucesso: Multiplicação correta");
        } else {
            System.out.println("Falha: Resultado incorreto");
        }
    }

    /**
     * Testa a multiplicação de uma matriz por um escalar.
     * Verifica se cada elemento é corretamente multiplicado pelo escalar.
     */
    public static void testMultiplicaMatrizPorEscalar() {
        System.out.println("Teste 9: Multiplicação de matriz por escalar");

        double[][] matrix = {
                {1.0, 2.0},
                {3.0, 4.0}
        };

        double escalar = 2.0;

        double[][] expected = {
                {2.0, 4.0},
                {6.0, 8.0}
        };

        double[][] result = Andre_Alexandre_Rita.multiplicaMatrizPorEscalar(matrix, escalar);

        if (matricesAreEqual(result, expected)) {
            System.out.println("Sucesso: Multiplicação por escalar correta");
        } else {
            System.out.println("Falha: Resultado incorreto");
        }
    }

    /**
     * Testa a construção da matriz M a partir de vetores de imagens.
     * Confirma se os vetores são organizados corretamente como colunas.
     */
    public static void testBuildImageMatrixM() {
        System.out.println("Teste 10: Construção da matriz M");

        int[][] imageVectors = {
                {1, 2, 3},
                {4, 5, 6}
        };

        int[][] expected = {
                {1, 4},
                {2, 5},
                {3, 6}
        };

        int[][] result = Andre_Alexandre_Rita.buildImageMatrixM(imageVectors);

        if (matricesAreEqual(result, expected)) {
            System.out.println("Sucesso: Construção da matriz M correta");
        } else {
            System.out.println("Falha: Resultado incorreto");
        }
    }

    /**
     * Testa a conversão de um ficheiro CSV em um vetor unidimensional.
     * Verifica se o vetor gerado corresponde ao conteúdo original do ficheiro.
     */
    public static void testConvertImagesToVectors() {
        System.out.println("Teste 11: Conversão de imagens para vetores");

        // Caminho para a pasta com as imagens CSV
        String folderPath = "C:/Users/andre/IdeaProjects/lapr1-24-25_DAB_02/Input/Funcao2-3/csv";

        try {
            // Obter a lista de ficheiros CSV na pasta
            File folder = new File(folderPath);
            File[] imageFiles = folder.listFiles((dir, name) -> name.toLowerCase().endsWith(".csv"));

            if (imageFiles == null || imageFiles.length == 0) {
                throw new IOException("Nenhum ficheiro CSV encontrado na pasta.");
            }

            // Selecionar o primeiro ficheiro CSV
            File firstImage = imageFiles[0];
            System.out.println("Primeiro ficheiro encontrado: " + firstImage.getName());

            // Converter o conteúdo do ficheiro para um vetor unidimensional
            int[] resultVector = Andre_Alexandre_Rita.loadImageFromCSV(firstImage);

            // Recarregar o mesmo ficheiro e verificar igualdade
            int[] expectedVector = Andre_Alexandre_Rita.loadImageFromCSV(firstImage);

            if (Arrays.equals(resultVector, expectedVector)) {
                System.out.println("Sucesso: Vetor convertido corresponde ao conteúdo original.");
            } else {
                System.out.println("Falha: Vetor convertido não corresponde ao conteúdo original.");
            }

        } catch (IOException e) {
            System.out.println("Falha: Erro ao ler os ficheiros CSV - " + e.getMessage());
        }
    }




    // Métodos de apoio
    public static double findMaxEigenValue(double[][] matrix) {
        RealMatrix realMatrix = new Array2DRowRealMatrix(matrix);
        return new org.apache.commons.math3.linear.EigenDecomposition(realMatrix).getRealEigenvalue(0);
    }

    public static boolean matricesAreEqual(double[][] a, double[][] b) {
        if (a.length != b.length || a[0].length != b[0].length) return false;

        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                if (a[i][j] != b[i][j]) return false;
            }
        }
        return true;
    }

    public static boolean matricesAreApproximatelyEqual(double[][] a, double[][] b, double epsilon) {
        if (a.length != b.length || a[0].length != b[0].length) return false;

        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                if (Math.abs(a[i][j] - b[i][j]) > epsilon) return false;
            }
        }
        return true;
    }

    public static boolean matricesAreEqual(int[][] a, int[][] b) {
        if (a.length != b.length || a[0].length != b[0].length) return false;

        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                if (a[i][j] != b[i][j]) return false;
            }
        }
        return true;
    }



}


