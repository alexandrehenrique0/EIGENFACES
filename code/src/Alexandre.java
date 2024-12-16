import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Arrays;public class Alexandre {
    // Calcula o Erro Absoluto Médio (EAM)
    public static void main(String[] args) {
        // Matriz de exemplo (3x3)
        double[][] Z = {
                {2.0, 0.0, 0.0},
                {0.0, 3.0, 0.0},
                {0.0, 0.0, 4.0}
        };

        int N = Z[0].length;
        int k = 2;

        double[] vetorMedio = vetorMedio(Z);

        System.out.println("Coluna Media");
        for (int i = 0; i < vetorMedio.length; i++) {
            System.out.println(vetorMedio[i]);

        }

        double[][] desviosA = matrixDesvios(Z,vetorMedio);
        double[][] desviosAT = transpostaMatriz(desviosA);

        print_Matrix(desviosA, " desvioA");
        print_Matrix(desviosAT, " desvioAT");

        double[][] covariancia = covariancias(desviosA,N);
        print_Matrix(covariancia, "matriz C : covariancia ");

        double[][] valProATxA = valoresPropriosATxA(desviosA,desviosAT);
        print_Matrix(valProATxA, "VALORES A^t . A");

        double[][] vi = vetoresPropriosATxA(desviosA,desviosAT);
        print_Matrix(vi, "Vetores A^t . A");

        double[][] vetProAxAT = vetoresPropriosAxAT(desviosA, vi);
        print_Matrix(vetProAxAT, "Vetores A . vi");

        double[][] valProC = valoresPropriosC(valProATxA,N);
        print_Matrix(valProC, "VALORES C");

        double[][] vetProC = vetProAxAT;
        print_Matrix(vetProC, "Vetores C");

        double[][] vetNormalizados = normalizarVetores(vetProC);
        print_Matrix(vetNormalizados, "Vetores Normalizados");

        double[] imagemNova = {2.0, 0.0, 0.0};

        double[][] phi = phiNova(imagemNova, vetorMedio);
        for (int i = 0; i < phi.length; i++) {
            System.out.println(phi[i][0]);
        };

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

    private static void print_Line(int length, String pattern) {
        for (int i = 0; i < length; i++) {
            System.out.print(pattern);
        }
        System.out.println();
    }

    private static void print_Matrix(double[][] matrixToPrint, String matrixName) {
        System.out.println("Matriz: " + matrixName + " ↓");
        print_Line(matrixToPrint[0].length, "");

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
    public static double[] vetorMedio(double[][] matrix) {
        int linhas = matrix.length, colunas = matrix[0].length;
        double[] media = new double[linhas];

        for (int i = 0; i < linhas; i++) {
            for (int j = 0; j < colunas; j++) {
                media[i] += matrix[i][j];
            }
            media[i] /= colunas;
        }
        return media;
    }
    public static boolean testColunaMedia(double[] matrix, double[] expectedColMedia) {
        return Arrays.equals(matrix, expectedColMedia);
    }

    public static double[][] matrixDesvios(double[][] matrix, double[] media) {
        int linhas = matrix.length, colunas = matrix[0].length;
        double[][] centralizada = new double[linhas][colunas];

        for (int i = 0; i < linhas; i++) {
            for (int j = 0; j < colunas; j++) {
                centralizada[i][j] = matrix[i][j] - media[i];
            }
        }
        return centralizada;

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

    public static double[][] vetoresPropriosAxAT(double[][] A,double[][] vi) {
        double[][] AxVI = multiplicaMatrizes(A,vi);
        return AxVI;
    }

    public static double[][] valoresPropriosAxAT(double[][] A,double[][] AT) {
        double[][] AxAT = multiplicaMatrizes(A,AT);
        EigenDecomposition eigenDecomposition = decomposeMatrix(AxAT);
        RealMatrix D = eigenDecomposition.getD();
        double[][] dArray = D.getData();
        adjustPrecision(dArray, 0.001);
        return dArray;
    }

    //ignora valores negativos muiiiito baixos, para nao obter -0.0
    private static double[][] adjustPrecision(double[][] matrix, double ignoreMinorValues) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (Math.abs(matrix[i][j]) < ignoreMinorValues) {
                    matrix[i][j] = 0.0;
                }
            }
        }
        return matrix;
    }

    public static double[][] valoresPropriosC(double[][] valoresPropriosAxAT, int N) {
        double[][] lambdai = multiplicaMatrizPorEscalar(valoresPropriosAxAT,1.0/N);
        return lambdai;
    }

    //------------------------------------------------------
    //5-----------------------------------------------------
    public static double[][] normalizarVetores(double[][] vetoresPropriosC) {
        for (int i = 0; i < vetoresPropriosC[0].length; i++) {
            double[] vetorAtual = new double[vetoresPropriosC.length];

            for (int j = 0; j < vetoresPropriosC.length; j++) {
                vetorAtual[j] = vetoresPropriosC[j][i];
            }

            double norm = calculateNorm(vetorAtual);

            for (int j = 0; j < vetoresPropriosC.length; j++) {
                vetoresPropriosC[j][i] = vetorAtual[j] / norm;
            }

        }
        return vetoresPropriosC;
    }


    private static double calculateNorm(double[] vector) {
        double sum = 0;
        for (int i = 0; i < vector.length; i++) {
            sum = sum + vector[i] * vector[i];
        }
        return Math.sqrt(sum);
    }

    //------------------------------------------------------
    //6-----------------------------------------------------
    public static double[][] phiNova(double[] imagemNova, double[] media) {
        int dimensao = imagemNova.length;
        double[][] phi = new double[dimensao][1];

        for (int i = 0; i < dimensao; i++) {
            phi[i][0] = imagemNova[i] - media[i];
        }

        return phi;
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