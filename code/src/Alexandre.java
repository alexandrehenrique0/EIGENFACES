import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

public class Alexandre {
    // Calcula o Erro Absoluto Médio (EAM)
    public static void main(String[] args) {
        // Matriz de exemplo (3x3)
        double[][] A = {
                {1.0, 2.0, 3.0},
                {4.0, 5.0, 6.0},
                {7.0, 8.0, 9.0}
        };

        double[][] colunaMedia = colunaMedia(A);
        System.out.println("Média das colunas:");
        for (int i = 0; i < colunaMedia.length; i++) {
            System.out.println(colunaMedia[i][0]);
        }
        double[][] desvios = calculoDesvios(A,colunaMedia);
        for (int i = 0; i < desvios.length; i++) {
            for (int j = 0; j < desvios[i].length; j++) {
                System.out.print(desvios[i][j] + " ");
            }
            System.out.println();
        }
        RealMatrix covariancia = covariancias(desvios);
        System.out.println(covariancia);

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

    public static boolean testCalculateEAM(double[][] A, double[][] Ak, double expectedEAM) {
        // Chama a função que calcula o EAM
        double calculatedEAM = calculateEAM(A, Ak);

        // Compara o valor calculado com o esperado
        return calculatedEAM == expectedEAM;
    }

    public static double[][] colunaMedia(double[][] matrix) {
        int linhas = matrix.length;
        int colunas = matrix[0].length;
        double[][] mediaColuna = new double[linhas][1];
        for (int i = 0; i < linhas; i++) {
            double soma = 0;
            for (int j = 0; j < colunas; j++) {
                soma = soma + matrix[i][j]; // Soma os elementos de cada linha
            }
            mediaColuna[i][0] = soma/colunas; // Calcula a média e armazena na matriz coluna
        }
        return mediaColuna;
    }
    public static boolean testColunaMedia(double[] matrix, double[] expectedColMedia) {
        return matrix == expectedColMedia;
    }

    public static double[][] calculoDesvios(double[][] matrix, double[][] colunaMedia) {
        int colunas = matrix[0].length;
        int linhas = matrix.length;
        double[][] desvios = new double[linhas][colunas];

        for (int j = 0; j < colunas; j++) {
            for (int i = 0; i < linhas; i++) {
                desvios[i][j] = matrix[i][j] - colunaMedia[i][0]; // Calcula o desvio de cada elemento da matriz em relação à média da linha correspondente
            }

        }
        return desvios;
    }

    public static boolean testCalDesvios(double[][] matrix, double[][] expectedDesvios) {
        return matrix == expectedDesvios;
    }

    public static RealMatrix covariancias(double[][] matrix) {
        int N = matrix.length;
        Array2DRowRealMatrix A = new Array2DRowRealMatrix(matrix);
        RealMatrix AT = A.transpose();
        RealMatrix C = A.multiply(AT).scalarMultiply(1.0/N);
        return C;
    }

    public static boolean testCovariancia(RealMatrix C, RealMatrix expectedC) {
        return C == expectedC;
    }
}