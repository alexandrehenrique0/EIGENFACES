import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import java.util.Arrays;
import java.util.Random;

public class Andre {

    /**
     * Método para normalizar os vetores próprios de uma matriz.
     * @param eigenVectors Matriz cujas colunas representam os vetores próprios.
     *                     A matriz é alterada diretamente com os vetores normalizados.
     */
    public static void normalizeEigenVectors(RealMatrix eigenVectors) {
        double[][] data = eigenVectors.getData(); // Extrair os dados como array bidimensional

        for (int i = 0; i < data[0].length; i++) { // Iterar por cada coluna
            double norm = 0.0;

            // Calcular a norma da coluna
            for (int j = 0; j < data.length; j++) {
                norm += data[j][i] * data[j][i];
            }
            norm = Math.sqrt(norm);

            // Verificar se a norma é válida antes de normalizar
            if (norm > 0) {
                for (int j = 0; j < data.length; j++) {
                    data[j][i] /= norm; // Normalizar cada valor
                }
            } else {
                // Reportar que a coluna é nula e deixá-la inalterada
                System.out.println("Vetor Próprio " + (i + 1) + " possui apenas zeros. Não foi normalizado.");
            }

            // Exibir o vetor normalizado ou inalterado
            System.out.print("Vetor Próprio Normalizado " + (i + 1) + ": [");
            for (int j = 0; j < data.length; j++) {
                System.out.print(data[j][i] + (j < data.length - 1 ? ", " : ""));
            }
            System.out.println("]");
        }

        // Substituir os dados normalizados na matriz original
        eigenVectors.setSubMatrix(data, 0, 0);
    }




    /**
     * Método para calcular a norma (L2) de um vetor.
     *
     * @param vector Vetor a ser normalizado
     * @return Norma do vetor
     */
    private static double calculateNorm(double[] vector) {
        double sum = 0;
        for (double v : vector) {
            sum += v * v;
        }
        return Math.sqrt(sum);
    }

    /**
     * Método para imprimir um vetor na consola.
     * @param vector Vetor a ser impresso
     */
    private static void printVector(double[] vector) {
        for (double v : vector) {
            System.out.print(v + " ");
        }
        System.out.println();
    }

    /**
     * Método para calcular a matriz de covariância a partir de uma matriz A.
     * @param A Matriz de entrada
     * @return Matriz de covariância
     */
    public static RealMatrix calculateCovarianceMatrix(RealMatrix A) {
        int numRows = A.getRowDimension();
        RealMatrix meanMatrix = calculateMeanMatrix(A);

        // Subtrair a média de cada linha da matriz A
        RealMatrix centeredMatrix = A.subtract(meanMatrix);

        // Calcular C = (1 / N) * A * A^T
        RealMatrix covarianceMatrix = centeredMatrix.multiply(centeredMatrix.transpose()).scalarMultiply(1.0 / numRows);
        return covarianceMatrix;
    }

    /**
     * Método para calcular a matriz de média (µ) de uma matriz.
     * @param A Matriz de entrada
     * @return Matriz contendo o vetor média repetido
     */
    private static RealMatrix calculateMeanMatrix(RealMatrix A) {
        int numRows = A.getRowDimension();
        int numCols = A.getColumnDimension();
        double[] mean = new double[numCols];

        // Calcular a média de cada coluna
        for (int j = 0; j < numCols; j++) {
            double sum = 0;
            for (int i = 0; i < numRows; i++) {
                sum += A.getEntry(i, j);
            }
            mean[j] = sum / numRows;
        }

        // Criar a matriz de média repetida
        double[][] meanMatrixData = new double[numRows][numCols];
        for (int i = 0; i < numRows; i++) {
            System.arraycopy(mean, 0, meanMatrixData[i], 0, numCols);
        }

        return new Array2DRowRealMatrix(meanMatrixData);
    }

    /**
     * Método principal para testes.
     */
    public static void testMain() {
        // Gerar uma matriz aleatória simétrica com tamanho até 256x256
        int size = 5; // Substituir pelo tamanho desejado
        RealMatrix testMatrix = generateRandomSymmetricMatrix(size);

        System.out.println("Testing covariance matrix calculation and normalization with random symmetric matrix of size " + size + ":");

        // Exibir a matriz gerada
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                System.out.print(testMatrix.getEntry(i, j) + " ");
            }
            System.out.println();
        }

        // Calcular a matriz de covariância (podem apagar se não estiver em conformidade, foi só para testes,implementem o do Alexandre)
        RealMatrix covarianceMatrix = calculateCovarianceMatrix(testMatrix);
        System.out.println("Covariance matrix:");
        for (int i = 0; i < covarianceMatrix.getRowDimension(); i++) {
            for (int j = 0; j < covarianceMatrix.getColumnDimension(); j++) {
                System.out.print(covarianceMatrix.getEntry(i, j) + " ");
            }
            System.out.println();
        }

        // Normalizar os eigenfaces da matriz de covariância
        EigenDecomposition decomposition = new EigenDecomposition(covarianceMatrix);
        normalizeEigenVectors(decomposition.getV());
    }

    /**
     * Método para gerar uma matriz simétrica aleatória com tamanho até 256x256.
     * @param size Tamanho da matriz (deve ser no máximo 256)
     * @return Matriz simétrica gerada
     */
    private static RealMatrix generateRandomSymmetricMatrix(int size) {
        if (size > 256) {
            throw new IllegalArgumentException("Matrix size must be 256 or smaller.");
        }

        double[][] matrixData = new double[size][size];
        Random random = new Random();

        for (int i = 0; i < size; i++) {
            for (int j = i; j < size; j++) {
                double value = random.nextDouble() * 100; // Valores aleatórios entre 0 e 100
                matrixData[i][j] = value;
                matrixData[j][i] = value; // Garantir simetria
            }
        }

        return new Array2DRowRealMatrix(matrixData);
    }
}
 /** Para testar na main apenas meter
  Andre.testMain();
  */