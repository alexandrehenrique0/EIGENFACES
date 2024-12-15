import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;


/**
 * Classe que integra os métodos desenvolvidos por André, Alexandre e Rita.
 * Esta classe unifica todas as funcionalidades necessárias para processamento
 * de Eigenfaces, como leitura de imagens, cálculo de covariância e normalização.
 *
 * Integração e unificação dos métodos realizada por: André.
 */

public class Andre_Alexandre_Rita {
    public static void main(String[] args) {
        // Caminho da pasta com os ficheiros CSV de entrada
        String csvPath = "Input/Funcao2-3/csv";

        // Método principal que integra todas as funcionalidades
        try {
            integrateEigenfaceProcesses(csvPath);
        } catch (IOException e) {
            System.err.println("Erro ao processar as imagens: " + e.getMessage());
        }
    }

    /**
     * Método que integra os processos dos colegas: Alexandre, Rita e André
     * @param csvPath Caminho para os ficheiros CSV das imagens
     */
    public static void integrateEigenfaceProcesses(String csvPath) throws IOException {
        // Leitura e conversão das imagens (Método da Rita)
        int[][] imageVectors = convertImagesToVectors(csvPath);
        int[][] imageMatrixM = buildImageMatrixM(imageVectors);
        System.out.println("=== Matriz de Imagens M ===");
        printMatrix(convertIntToDouble(imageMatrixM));

        // Transformação para double[][]
        double[][] imageMatrix = convertIntToDouble(imageMatrixM);

        // Cálculo da Coluna Média (Alexandre)
        double[][] colunaMedia = colunaMedia(imageMatrix);
        System.out.println("=== Coluna Média ===");
        printMatrix(colunaMedia);

        // Centralização das imagens
        double[][] desviosA = matrixDesvios(imageMatrix, colunaMedia);
        System.out.println("=== Matriz de Desvios (A) ===");
        printMatrix(desviosA);

        // Cálculo da Covariância (Alexandre)
        double[][] covariancia = covariancias(desviosA, imageMatrix[0].length);
        System.out.println("=== Matriz de Covariância ===");
        printMatrix(covariancia);

        // Cálculo dos valores e vetores próprios (André)
        RealMatrix covarianceMatrix = new Array2DRowRealMatrix(covariancia);
        EigenDecomposition eigenDecomposition = new EigenDecomposition(covarianceMatrix);
        RealMatrix eigenVectors = eigenDecomposition.getV();

        // Normalização dos vetores próprios (André)
        System.out.println("=== Vetores Próprios Normalizados ===");
        normalizeEigenVectors(eigenVectors);
    }

    /**
     * Converte uma matriz de int[][] para double[][]
     * Método desenvolvido por: Rita
     * @param matrix Matriz de entrada no formato int[][]
     * @return Matriz convertida para double[][]
     */
    public static double[][] convertIntToDouble(int[][] matrix) {
        double[][] result = new double[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                result[i][j] = matrix[i][j];
            }
        }
        return result;
    }

    /**
     * Calcula a coluna média de uma matriz
     * Método desenvolvido por: Alexandre
     * @param matrix Matriz de entrada
     * @return Uma Matriz que representa a média das colunas
     */
    public static double[][] colunaMedia(double[][] matrix) {
        int linhas = matrix.length;
        int colunas = matrix[0].length;
        double[][] mediaColuna = new double[linhas][1];
        for (int i = 0; i < linhas; i++) {
            double soma = 0;
            for (int j = 0; j < colunas; j++) {
                soma += matrix[i][j];
            }
            mediaColuna[i][0] = soma / colunas;
        }
        return mediaColuna;
    }

    /**
     * Centraliza a matriz em relação à coluna média
     * Método desenvolvido por: Alexandre
     * @param matrix Matriz original
     * @param colunaMedia Coluna média calculada
     * @return Matriz centralizada
     */
    public static double[][] matrixDesvios(double[][] matrix, double[][] colunaMedia) {
        int linhas = matrix.length;
        int colunas = matrix[0].length;
        double[][] desvios = new double[linhas][colunas];
        for (int i = 0; i < linhas; i++) {
            for (int j = 0; j < colunas; j++) {
                desvios[i][j] = matrix[i][j] - colunaMedia[i][0];
            }
        }
        return desvios;
    }

    /**
     * Calcula a matriz de covariância
     * Método desenvolvido por: Alexandre
     * @param A Matriz centralizada
     * @param N Quantidade de imagens
     * @return Matriz de covariância
     */
    public static double[][] covariancias(double[][] A, int N) {
        double[][] AT = transpostaMatriz(A);
        double[][] AAT = multiplicaMatrizes(A, AT);
        return multiplicaMatrizPorEscalar(AAT, 1.0 / N);
    }

    /**
     * Normaliza os vetores próprios de uma matriz.
     * Cada coluna da matriz 'eigenVectors' é tratada como um vetor próprio
     * e será normalizada (dividida pela sua norma Euclidiana).
     *
     * Caso a norma de um vetor seja zero, ele é identificado como um vetor nulo
     * e permanece inalterado.
     *
     * Método desenvolvido por: André
     *
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
     * Calcula a norma de um vetor
     * Método desenvolvido por: André
     * @param vector Vetor de entrada
     * @return Norma calculada
     */
    private static double calculateNorm(double[] vector) {
        double sum = 0;
        for (double v : vector) {
            sum += v * v;
        }
        return Math.sqrt(sum);
    }

    /**
     * Lê e converte imagens CSV para vetores
     * Método desenvolvido por: Rita
     * @param imageFolderLocation Caminho da pasta contendo os ficheiros CSV
     * @return Matriz de vetores de imagens
     * @throws IOException Se ocorrer erro na leitura dos ficheiros
     */
    public static int[][] convertImagesToVectors(String imageFolderLocation) throws IOException {
        File folder = new File(imageFolderLocation);
        File[] imageFiles = folder.listFiles((dir, name) -> name.toLowerCase().endsWith(".csv"));

        if (imageFiles == null || imageFiles.length == 0) {
            throw new IOException("Nenhuma imagem CSV encontrada na pasta.");
        }

        int[][] imageVectors = new int[imageFiles.length][];
        for (int i = 0; i < imageFiles.length; i++) {
            imageVectors[i] = loadImageFromCSV(imageFiles[i]);
        }
        return imageVectors;
    }

    /**
     * Lê o conteúdo de um ficheiro CSV e converte-o num vetor de inteiros.
     * Método desenvolvido por: Rita
     * @param imageFile Ficheiro CSV de entrada
     * @return Vetor contendo os dados do ficheiro
     * @throws IOException Se ocorrer um erro na leitura do ficheiro
     */
    public static int[] loadImageFromCSV(File imageFile) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(imageFile));
        String line;
        int numRows = 0;
        int numCols = 0;

        // Contar o número de linhas e colunas
        while ((line = reader.readLine()) != null) {
            numCols = line.split(",").length;
            numRows++;
        }
        reader.close();

        int[] imageVector = new int[numRows * numCols];
        reader = new BufferedReader(new FileReader(imageFile));
        int index = 0;

        // Ler os dados do ficheiro CSV e armazenar no vetor
        while ((line = reader.readLine()) != null) {
            for (String pixel : line.split(",")) {
                imageVector[index++] = Integer.parseInt(pixel.trim());
            }
        }
        reader.close();
        return imageVector;
    }


    /**
     * Constrói uma matriz M a partir dos vetores de imagens
     * Método desenvolvido por: Rita
     * @param imageVectors Vetores de imagens
     * @return Matriz M
     */
    public static int[][] buildImageMatrixM(int[][] imageVectors) {
        int numPixels = imageVectors[0].length;
        int numImages = imageVectors.length;
        int[][] imageMatrixM = new int[numPixels][numImages];
        for (int i = 0; i < numImages; i++) {
            for (int j = 0; j < numPixels; j++) {
                imageMatrixM[j][i] = imageVectors[i][j];
            }
        }
        return imageMatrixM;
    }


/**
     * Imprime uma matriz no formato legível
     * @param matrix Matriz a ser impressa
     */
    public static void printMatrix(double[][] matrix) {
        for (double[] row : matrix) {
            System.out.println(Arrays.toString(row));
        }
    }

    /**
     * Transpõe uma matriz
     * @param A Matriz original
     * @return Matriz transposta
     */
    public static double[][] transpostaMatriz(double[][] A) {
        int linhas = A.length;
        int colunas = A[0].length;
        double[][] transposta = new double[colunas][linhas];
        for (int i = 0; i < linhas; i++) {
            for (int j = 0; j < colunas; j++) {
                transposta[j][i] = A[i][j];
            }
        }
        return transposta;
    }

    /**
     * Multiplica duas matrizes
     * @param A Matriz A
     * @param B Matriz B
     * @return Resultado da multiplicação
     */
    public static double[][] multiplicaMatrizes(double[][] A, double[][] B) {
        int linhas = A.length;
        int colunas = B[0].length;
        double[][] resultado = new double[linhas][colunas];
        for (int i = 0; i < linhas; i++) {
            for (int j = 0; j < colunas; j++) {
                for (int k = 0; k < A[0].length; k++) {
                    resultado[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return resultado;
    }

    /**
     * Multiplica uma matriz por um escalar
     * @param A Matriz original
     * @param escalar Valor escalar
     * @return Matriz resultante
     */
    public static double[][] multiplicaMatrizPorEscalar(double[][] A, double escalar) {
        int linhas = A.length;
        int colunas = A[0].length;
        double[][] resultado = new double[linhas][colunas];
        for (int i = 0; i < linhas; i++) {
            for (int j = 0; j < colunas; j++) {
                resultado[i][j] = A[i][j] * escalar;
            }
        }
        return resultado;
    }

}
