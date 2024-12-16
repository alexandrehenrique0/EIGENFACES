import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Arrays;

public class Alexandre {

    /**
     * Constrói a matriz M onde cada coluna representa uma imagem.
     */
    public static double[][] buildImageMatrixM ( int[][] imageVectors){
        int numPixels = imageVectors[0].length;
        int numImages = imageVectors.length;
        double[][] matrixM = new double[numPixels][numImages];

        for (int i = 0; i < numImages; i++) {
            for (int j = 0; j < numPixels; j++) {
                matrixM[j][i] = imageVectors[i][j];
            }
        }
        return matrixM;
    }

    /**
     * Calcula o vetor médio μ.
     */
    public static double[] calcularVetorMedio ( double[][] matriz){
        int linhas = matriz.length, colunas = matriz[0].length;
        double[] media = new double[linhas];

        for (int i = 0; i < linhas; i++) {
            for (int j = 0; j < colunas; j++) {
                media[i] += matriz[i][j];
            }
            media[i] /= colunas;
        }
        return media;
    }

    /**
     * Centraliza a matriz em relação ao vetor médio.
     */
    public static double[][] centralizarMatriz ( double[][] matriz, double[] media){
        int linhas = matriz.length, colunas = matriz[0].length;
        double[][] centralizada = new double[linhas][colunas];

        for (int i = 0; i < linhas; i++) {
            for (int j = 0; j < colunas; j++) {
                centralizada[i][j] = matriz[i][j] - media[i];
            }
        }
        return centralizada;
    }

    /**
     * Calcula a matriz de covariância.
     */
    public static double[][] calcularCovariancia ( double[][] matrizA){
        int linhas = matrizA.length, colunas = matrizA[0].length;
        double[][] covariancia = new double[linhas][linhas];

        for (int i = 0; i < linhas; i++) {
            for (int j = 0; j < linhas; j++) {
                double soma = 0;
                for (int k = 0; k < colunas; k++) {
                    soma += matrizA[i][k] * matrizA[j][k];
                }
                covariancia[i][j] = soma / colunas;
            }
        }
        return covariancia;
    }

}
/*    // Calcula o Erro Absoluto Médio (EAM)
    public static void main(String[] args) {
        // Matriz de exemplo (3x3)
        double[][] Z = {
                {1.0, 2.0, 3.0},
                {2.0, 1.0, 6.0},
                {3.0, 6.0, 1.0}
        };

        int N = Z[0].length;
        int k = 2;

        double[][] colunaMedia = colunaMedia(Z);

        printMatrix(colunaMedia, "Coluna Média");

        double[][] desviosA = matrixDesvios(Z,colunaMedia);
        double[][] desviosAT = transpostaMatriz(desviosA);

        printMatrix(desviosA, " desvioA");
        printMatrix(desviosAT, " desvioAT");

        printMatrix(desviosA, "Matriz Z = desvioA");

        double[][] covariancia = covariancias(desviosA,N);
        printMatrix(covariancia, "matriz C : covariancia ");

        double[][] valProATxA = valoresPropriosATxA(desviosA,desviosAT);
        printMatrix(valProATxA, "VALORES A^t . A");

        double[][] vetProATxA = vetoresPropriosATxA(desviosA,desviosAT);
        printMatrix(vetProATxA, "Vetores A^t . A");

        double[][] vetProAxAT = vetoresPropriosAxAT(desviosA,desviosAT);
        printMatrix(vetProAxAT, "Vetores A . A^t");


        // algo esta mal nessa funcao
        double[][] valProAxAt = valoresPropriosAxAT(desviosA,desviosAT);
        printMatrix(valProAxAt, "VALORES A^t . A");

        double[][] valProC = valoresPropriosC(valProATxA,N);
        printMatrix(valProC, "VALORES C");

        double[][] vetProC = vetProAxAT;
        printMatrix(vetProC, "Vetores C");

        double[][] vetNormalizados = normalizarVetores(Z);
        printMatrix(vetNormalizados, "Vetores Normalizados");

    }


    /**
     * Calcula o Erro Absoluto Médio (EAM) entre duas matrizes.
     * @param matrizOriginal Matriz original.
     * @param matrizReconstruida Matriz reconstruída.
     * @return O valor do erro absoluto médio entre as duas matrizes.
     * @throws IllegalArgumentException Se as dimensões das matrizes não coincidirem.
     */
   /* public static double calcularErroAbsolutoMedio(double[][] matrizOriginal, double[][] matrizReconstruida) {
        int linhas = matrizOriginal.length;
        int colunas = matrizOriginal[0].length;

        // Verifica se as dimensões das matrizes são compatíveis
        if (linhas != matrizReconstruida.length || colunas != matrizReconstruida[0].length) {
            throw new IllegalArgumentException("As dimensões das matrizes não coincidem.");
        }

        double somaErroAbsoluto = 0;

        // Percorre cada elemento da matriz
        for (int i = 0; i < linhas; i++) {
            for (int j = 0; j < colunas; j++) {
                somaErroAbsoluto += Math.abs(matrizOriginal[i][j] - matrizReconstruida[i][j]);
            }
        }

        // Calcula o erro médio
        return somaErroAbsoluto / (linhas * colunas);
    }



    private static void print_Line(int length, String pattern) {
        for (int i = 0; i < length; i++) {
            System.out.print(pattern);
        }
        System.out.println();
    }

    /**
     * Imprime uma matriz no formato legível.
     * @param matrix Matriz a ser impressa.
     */
 /*   public static void printMatrix(double[][] matrix) {
        if (matrix == null || matrix.length == 0) {
            System.out.println("A matriz está vazia ou é nula.");
            return;
        }

        for (double[] row : matrix) {
            System.out.println(Arrays.toString(row));
        }
    }

    /*public static boolean testCalculateEAM(double[][] A, double[][] Ak, double expectedEAM) {
        // Chama a função que calcula o EAM
       // double calculatedEAM = calculateEAM(A, Ak);

        // Compara o valor calculado com o esperado
      //  return calculatedEAM == expectedEAM;
 //   }

    /**
     * Calcula o vetor médio μ, a média de todas as colunas (imagens).
     * @param matrizM Matriz M onde cada coluna é uma imagem vetorizada.
     * @return Vetor médio μ.
     */
    /*public static double[] calcularVetorMedio(double[][] matrizM) {
        int linhas = matrizM.length;
        int colunas = matrizM[0].length;
        double[] vetorMedio = new double[linhas];

        for (int i = 0; i < linhas; i++) {
            double soma = 0;
            for (int j = 0; j < colunas; j++) {
                soma += matrizM[i][j];
            }
            vetorMedio[i] = soma / colunas;
        }
        return vetorMedio;
    }

    /**
     * Centraliza a matriz em relação ao vetor médio μ.
     * @param matrizM Matriz original M.
     * @param vetorMedio Vetor médio μ.
     * @return Matriz centralizada A.
     */
   /* public static double[][] centralizarMatriz(double[][] matrizM, double[] vetorMedio) {
        int linhas = matrizM.length;
        int colunas = matrizM[0].length;
        double[][] matrizA = new double[linhas][colunas];

        for (int j = 0; j < colunas; j++) { // Cada coluna (imagem)
            for (int i = 0; i < linhas; i++) { // Cada linha (pixel)
                matrizA[i][j] = matrizM[i][j] - vetorMedio[i];
            }
        }
        return matrizA;
    }

    /**
     * Calcula a matriz de covariâncias C = (1/N) * A * A^T.
     * @param matrizA Matriz centralizada A.
     * @return Matriz de covariâncias C.
     */
  /*  public static double[][] calcularCovariancia(double[][] matrizA) {
        int linhas = matrizA.length;
        int colunas = matrizA[0].length;
        double[][] covariancia = new double[linhas][linhas];

        // Produto A * A^T
        for (int i = 0; i < linhas; i++) {
            for (int j = 0; j < linhas; j++) {
                double soma = 0;
                for (int k = 0; k < colunas; k++) {
                    soma += matrizA[i][k] * matrizA[j][k];
                }
                covariancia[i][j] = soma / colunas; // Divide por N
            }
        }
        return covariancia;
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

    public static EigenDecomposition decomposeMatrix(double[][] arrayParaDecompor) {
        Array2DRowRealMatrix matrix = new Array2DRowRealMatrix(arrayParaDecompor);
        EigenDecomposition eigenDecomposition = new EigenDecomposition(matrix);

//        RealMatrix eigenVectors = eigenDecomposition.getV();
//        RealMatrix eigenValues = eigenDecomposition.getD();
//        RealMatrix eigenVectorsTranspose = eigenDecomposition.getVT();

        return eigenDecomposition;
    }

    /**
     * Transpõe uma matriz (troca linhas por colunas).
     * @param matrizOriginal Matriz original a ser transposta.
     * @return Matriz transposta.
     * @throws IllegalArgumentException Se a matriz for nula ou irregular.
     */
   /* public static double[][] transpostaMatriz(double[][] matrizOriginal) {
        if (matrizOriginal == null || matrizOriginal.length == 0) {
            throw new IllegalArgumentException("A matriz original não pode ser nula ou vazia.");
        }

        int linhas = matrizOriginal.length;
        int colunas = matrizOriginal[0].length;

        // Verifica se a matriz é regular
        for (double[] linha : matrizOriginal) {
            if (linha.length != colunas) {
                throw new IllegalArgumentException("A matriz fornecida é irregular.");
            }
        }

        double[][] matrizTransposta = new double[colunas][linhas];

        for (int i = 0; i < linhas; i++) {
            for (int j = 0; j < colunas; j++) {
                matrizTransposta[j][i] = matrizOriginal[i][j];
            }
        }

        return matrizTransposta;
    }


    /**
     * Multiplica duas matrizes.
     * @param A Primeira matriz.
     * @param B Segunda matriz.
     * @return Resultado da multiplicação de A por B.
     * @throws IllegalArgumentException Se as dimensões das matrizes não permitirem a multiplicação.
     */
   /* public static double[][] multiplicaMatrizes(double[][] A, double[][] B) {
        int linhasA = A.length;
        int colunasA = A[0].length;
        int linhasB = B.length;
        int colunasB = B[0].length;

        // Verifica se as dimensões permitem a multiplicação
        if (colunasA != linhasB) {
            throw new IllegalArgumentException("O número de colunas da matriz A deve ser igual ao número de linhas da matriz B.");
        }

        double[][] resultado = new double[linhasA][colunasB];

        for (int i = 0; i < linhasA; i++) {
            for (int j = 0; j < colunasB; j++) {
                for (int k = 0; k < colunasA; k++) {
                    resultado[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return resultado;
    }

    /**
     * Multiplica uma matriz por um valor escalar.
     * @param matrizOriginal Matriz original.
     * @param escalar Valor escalar a ser multiplicado.
     * @return Nova matriz resultante da multiplicação.
     * @throws IllegalArgumentException Se a matriz for nula ou vazia.
     */
   /* public static double[][] multiplicaMatrizPorEscalar(double[][] matrizOriginal, double escalar) {
        if (matrizOriginal == null || matrizOriginal.length == 0 || matrizOriginal[0].length == 0) {
            throw new IllegalArgumentException("A matriz original não pode ser nula ou vazia.");
        }

        int linhas = matrizOriginal.length;
        int colunas = matrizOriginal[0].length;
        double[][] resultado = new double[linhas][colunas];

        for (int i = 0; i < linhas; i++) {
            for (int j = 0; j < colunas; j++) {
                resultado[i][j] = matrizOriginal[i][j] * escalar;
            }
        }

        return resultado;
    }
*/