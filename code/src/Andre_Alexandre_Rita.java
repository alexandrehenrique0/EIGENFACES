import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.*;
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
        String csvPath = "Input/Funcao2-3/csv";

        try {
            // 1. Carregar e converter imagens CSV para matriz M
            int[][] imageVectors = convertImagesToVectors(csvPath);
            double[][] matrizM = buildImageMatrixM(imageVectors);
            printMatrix(matrizM, "Matriz M (Imagens Vetorizadas)");
            System.out.println("✅ Matriz M (Imagens Vetorizadas) construída com sucesso!");

            double[] vetorMedio = calcularVetorMedio(matrizM);
            System.out.println("\n=== Vetor Médio (μ) ===");
            System.out.println(Arrays.toString(vetorMedio));
            System.out.println("✅ Vetor médio calculado com sucesso!");

            double[][] matrizA = centralizarMatriz(matrizM, vetorMedio);
            //printMatrix(matrizA, "Matriz A (Imagens Centralizadas)");
            System.out.println("✅ Matriz A (Imagens Centralizadas) construída com sucesso!");

            double[][] matrizC = calcularCovariancia(matrizA);
            //printMatrix(matrizC, "Matriz de Covariâncias C");
            System.out.println("✅ Matriz de covariâncias C calculada com sucesso!");

            // 5. Calcular os vetores próprios e normalizá-los
            RealMatrix covarianceMatrix = new Array2DRowRealMatrix(matrizC);
            EigenDecomposition eigenDecomposition = new EigenDecomposition(covarianceMatrix);
            RealMatrix eigenVectors = eigenDecomposition.getV();
            normalizeEigenVectors(eigenVectors);
            System.out.println("✅ Vetores próprios calculados e normalizados com sucesso!");

        } catch (IOException e) {
            System.err.println("Erro ao processar os ficheiros CSV: " + e.getMessage());
        }
    }

        /**
         * Lê ficheiros CSV e converte-os em vetores unidimensionais.
         */
        public static int[][] convertImagesToVectors (String imageFolderLocation) throws IOException {
            File folder = new File(imageFolderLocation);

            File[] imageFiles = folder.listFiles();

            if (imageFiles == null || imageFiles.length == 0) {
                throw new FileNotFoundException("Nenhuma imagem CSV encontrada na pasta.");
            }

            int csvFileCount = 0;
            for (int i = 0; i < imageFiles.length; i++) {
                File file = imageFiles[i];
                if (file.isFile() && file.getName().toLowerCase().endsWith(".csv")) {
                    csvFileCount++;
                }
            }

            if (csvFileCount == 0) {
                throw new FileNotFoundException("Nenhum arquivo CSV encontrado.");
            }

            File[] csvFiles = new File[csvFileCount];
            int index = 0;

            for (int i = 0; i < imageFiles.length; i++) {
                File file = imageFiles[i];
                if (file.isFile() && file.getName().toLowerCase().endsWith(".csv")) {
                    csvFiles[index++] = file;
                }
            }

            System.out.printf("Encontrados %d arquivos CSV.%n", csvFiles.length);

            int[][] imageVectors = new int[csvFiles.length][];

            for (int i = 0; i < csvFiles.length; i++) {
                imageVectors[i] = loadImageFromCSV(csvFiles[i]);
            }

            return imageVectors;
        }
        public static int[] loadImageFromCSV (File imageFile) throws IOException {

            BufferedReader reader = new BufferedReader(new FileReader(imageFile));
            String line;
            int numRows = 0;
            int numCols = 0;

            while ((line = reader.readLine()) != null) {
                String[] pixels = line.split(",");
                numCols = pixels.length;
                numRows++;
            }
            reader.close();

            String[][] imageData = new String[numRows][numCols];

            reader = new BufferedReader(new FileReader(imageFile));
            int rowIndex = 0;
            while ((line = reader.readLine()) != null) {
                String[] pixels = line.split(",");
                for (int colIndex = 0; colIndex < pixels.length; colIndex++) {
                    imageData[rowIndex][colIndex] = pixels[colIndex].trim();
                }
                rowIndex++;
            }
            reader.close();

            int[] imageVector = new int[numRows * numCols];
            int index = 0;

            for (int col = 0; col < numCols; col++) {
                for (int row = 0; row < numRows; row++) {
                    imageVector[index++] = Integer.parseInt(imageData[row][col]);
                }
            }

            return imageVector;
        }

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
            int N = matrizA[0].length;
            double[][] AT = transpostaMatriz(matrizA);
            double[][] AAT = multiplicaMatrizes(matrizA,AT);
            return multiplicaMatrizPorEscalar(AAT,1.0/N);
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

        /**
         * Normaliza os vetores próprios.
         */
        public static void normalizeEigenVectors (RealMatrix eigenVectors){
            double[][] data = eigenVectors.getData();

            for (int i = 0; i < data[0].length; i++) { // Iterar por cada coluna (vetor próprio)
                double norm = 0.0;

                // Calcular a norma Euclidiana do vetor (coluna)
                for (int j = 0; j < data.length; j++) {
                    norm += data[j][i] * data[j][i];
                }
                norm = Math.sqrt(norm);

                // Verificar se a norma é maior que zero antes de normalizar
                if (norm > 0) {
                    for (int j = 0; j < data.length; j++) {
                        data[j][i] /= norm; // Normalizar o vetor
                    }
                } else {
                    System.out.println("Vetor Próprio " + (i + 1) + " é nulo e não foi normalizado.");
                }
            }

            eigenVectors.setSubMatrix(data, 0, 0);
        }

        /**
         * Imprime uma matriz.
         */
        public static void printMatrix ( double[][] matrix, String title){
            System.out.println("\n=== " + title + " ===");
            for (double[] row : matrix) {
                for (double value : row) {
                    System.out.print(String.format("%.4f ", value));
                }
                System.out.println();
            }
            System.out.flush(); // Força a saída imediata
        }
    }
