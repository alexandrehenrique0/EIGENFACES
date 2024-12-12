import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Arrays;

public class Alexandre {
    // Calcula o Erro Absoluto Médio (EAM)
    public static void main(String[] args) {
        // Matriz de exemplo (3x3)
        double[][] Z = {
                {1.0, 2.0, 3.0},
                {2.0, 1.0, 6.0},
                {3.0, 6.0, 1.0}
        };

        int N = Z[0].length;

        //GABRIEL

        EigenDecomposition eigenDecomposition = decomposeMatrix(Z);
        double[][] eigenVectorsArray = getEigenVectors(eigenDecomposition);
        double[][] eigenValuesArray = getEigenValues(eigenDecomposition);
        double[][] eigenVectorsTransposeArray = getEigenVectorsTranspose(eigenDecomposition);


        //GABRIEL

        double[][] colunaMedia = colunaMedia(Z);

        print_Matrix(colunaMedia, "Coluna Média");

        double[][] desviosA = matrixDesvios(Z,colunaMedia);
        double[][] desviosAT = transpostaMatriz(desviosA);

        print_Matrix(desviosA, " desvioA");
        print_Matrix(desviosAT, " desvioAT");

        print_Matrix(desviosA, "Matriz Z (?)");

        double[][] covariancia = covariancias(desviosA,N);
        print_Matrix(covariancia, "matriz C : covariancia ");

        double[][] valProATxA = valoresPropriosATxA(desviosA,desviosAT);
        print_Matrix(valProATxA, "VALORES A^t . A");

        double[][] vetProATxA = vetoresPropriosATxA(desviosA,desviosAT);
        print_Matrix(vetProATxA, "Vetores A^t . A");

        double[][] vetProAxAT = vetoresPropriosAxAT(desviosA,desviosAT);
        print_Matrix(vetProAxAT, "Vetores A . A^t");

        double[][] valProC = valoresPropriosC(valProATxA,N);
        print_Matrix(valProC, "VALORES C");

        double[][] vetProC = vetProAxAT;
        print_Matrix(vetProC, "Vetores C");
    }

    //GABRIEL
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
    //GABRIEL

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

    public static boolean testCalculateEAM(double[][] A, double[][] Ak, double expectedEAM) {
        // Chama a função que calcula o EAM
        double calculatedEAM = calculateEAM(A, Ak);

        // Compara o valor calculado com o esperado
        return calculatedEAM == expectedEAM;
    }

    //3-----------------------------------------------------
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
        return Arrays.equals(matrix, expectedColMedia);
    }

    public static double[][] matrixDesvios(double[][] matrix, double[][] colunaMedia) {
        int colunas = matrix[0].length;
        int linhas = matrix.length;
        double[][] desvios = new double[linhas][colunas];

        for (int j = 0; j < colunas; j++) {
            for (int i = 0; i < linhas; i++) {
                desvios[i][j] =  colunaDesvio(matrix[i][j],colunaMedia[i][0]);
            }

        }
        return desvios;
    }

    public static double colunaDesvio(double valorMatrix, double valorColunaMedia) {

        double desvio = valorMatrix - valorColunaMedia; // Calcula o desvio de cada elemento da matriz em relação à média da linha correspondente

        return desvio;
    }

    /* public static double[][] calculoDesvios(double[][] matrix, double[][] colunaMedia) {
        int colunas = matrix[0].length;
        int linhas = matrix.length;
        double[][] desvios = new double[linhas][colunas];

        for (int j = 0; j < colunas; j++) {
            for (int i = 0; i < linhas; i++) {
                desvios[i][j] = matrix[i][j] - colunaMedia[i][0]; // Calcula o desvio de cada elemento da matriz em relação à média da linha correspondente
            }

        }
        return desvios;
    }*/

    public static boolean testCalDesvios(double[][] matrix, double[][] expectedDesvios) {
        return Arrays.equals(matrix, expectedDesvios);
    }

    public static double[][] covariancias(double[][] A,int N) {
        double[][] AT = transpostaMatriz(A);
        double[][] AAT = multiplicaMatrizes(A,AT);
        return multiplicaMatrizPorEscalar(AAT,1.0/N);
    }

    public static boolean testCovariancia(double[][] C, double[][] expectedC) {
        return C == expectedC;
    }



    //-------------------------------------------------------

    //4------------------------------------------------------
    public static double[][] valoresPropriosATxA(double[][] A,double[][] AT) {
        double[][] ATxA = multiplicaMatrizes(AT,A);
        EigenDecomposition eigenDecomposition = decomposeMatrix(ATxA);
        RealMatrix D = eigenDecomposition.getD();
        double[][] dArray = D.getData();
        return dArray;
    }

    public static double[][] vetoresPropriosATxA(double[][] A,double[][] AT) {
        double[][] ATxA = multiplicaMatrizes(AT,A);
        EigenDecomposition eigenDecomposition = decomposeMatrix(ATxA);
        RealMatrix V = eigenDecomposition.getV();
        double[][] vi = V.getData();
        return vi;
    }

    public static double[][] vetoresPropriosAxAT(double[][] A,double[][] AT) {
        double[][] AxAT = multiplicaMatrizes(A,AT);
        EigenDecomposition eigenDecomposition = decomposeMatrix(AxAT);
        RealMatrix V = eigenDecomposition.getV();
        double[][] ui = V.getData();
        return ui;
    }

    public static double[][] valoresPropriosC(double[][] valoresPropriosAxAT, int N) {
        double[][] lambdai = multiplicaMatrizPorEscalar(valoresPropriosAxAT,1.0/N);
        return lambdai;
    }

    //------------------------------------------------------

    public static EigenDecomposition decomposeMatrix(double[][] arrayParaDecompor) {
        Array2DRowRealMatrix matrix = new Array2DRowRealMatrix(arrayParaDecompor);
        EigenDecomposition eigenDecomposition = new EigenDecomposition(matrix);

//        RealMatrix eigenVectors = eigenDecomposition.getV();
//        RealMatrix eigenValues = eigenDecomposition.getD();
//        RealMatrix eigenVectorsTranspose = eigenDecomposition.getVT();

        return eigenDecomposition;
    }

    public static double[][] transpostaMatriz(double[][] matriz) {
        double[][] matrizTransposta = new double[matriz[0].length][matriz.length];
        for (int i = 0; i < matriz.length; i++) {
            for (int j = 0; j < matriz[0].length; j++) {
                matrizTransposta[j][i] = matriz[i][j];
            }
        }
        return matrizTransposta;
    }

    public static double[][] multiplicaMatrizes(double[][] matrizLeft, double[][] matrizRight) {
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

    public static double[][] multiplicaMatrizPorEscalar(double[][] matriz, double escalar) {
        double[][] matrizResultante = new double[matriz.length][matriz[0].length];
        for (int i = 0; i < matriz.length; i++) {
            for (int j = 0; j < matriz[0].length; j++) {
                matrizResultante[i][j] = matriz[i][j] * escalar;
            }
        }
        return matrizResultante;
    }
}