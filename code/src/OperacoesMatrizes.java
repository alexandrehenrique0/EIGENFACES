public class OperacoesMatrizes {
    
    public static double[][] somaMatrizes(double[][] matrizLeft, double[][] matrizRight) {
        double[][] matrizResultante = new double[matrizLeft.length][matrizLeft[0].length];
        for (int i = 0; i < matrizLeft.length; i++) {
            for (int j = 0; j < matrizLeft[0].length; j++) {
                matrizResultante[i][j] = matrizLeft[i][j] + matrizRight[i][j];
            }
        }
        return matrizResultante;
    }

    public static double[][] subtraiMatrizes(double[][] matrizLeft, double[][] matrizRight) {
        double[][] matrizResultante = new double[matrizLeft.length][matrizLeft[0].length];
        for (int i = 0; i < matrizLeft.length; i++) {
            for (int j = 0; j < matrizLeft[0].length; j++) {
                matrizResultante[i][j] = matrizLeft[i][j] - matrizRight[i][j];
            }
        }
        return matrizResultante;
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

    public static double[][] transpostaMatriz(double[][] matriz) {
        double[][] matrizTransposta = new double[matriz[0].length][matriz.length];
        for (int i = 0; i < matriz.length; i++) {
            for (int j = 0; j < matriz[0].length; j++) {
                matrizTransposta[j][i] = matriz[i][j];
            }
        }
        return matrizTransposta;
    }

}
