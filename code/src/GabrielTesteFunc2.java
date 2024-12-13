
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class GabrielTesteFunc2 {
    // Calcula o Erro Absoluto Médio (EAM)
    public static void main(String[] args) {
        String csvPath = "Input/TesteFuncao2-3/csv/image_001.csv";
        double[][] arrayOriginal = readCSVToArray(csvPath);

        // Dimensão e número de eigenfaces
        int N = arrayOriginal[0].length; // Número de colunas (64)
        int k = 64; // Número de eigenfaces (ajuste conforme necessário)

        // 1. Cálculo da média
        double[][] colunaMedia = colunaMedia(arrayOriginal);
        print_Matrix(colunaMedia, "Coluna Média");

        // 2. Centralização (Matriz de desvios)
        double[][] desviosA = matrixDesvios(arrayOriginal, colunaMedia);
        double[][] desviosAT = transpostaMatriz(desviosA);
        print_Matrix(desviosA, "Matriz de Desvios (A)");
        print_Matrix(desviosAT, "Matriz Transposta de Desvios (A^T)");

        // 3. Matriz de Covariância
        double[][] covariancia = covariancias(desviosA, N);
        print_Matrix(covariancia, "Matriz de Covariância (C)");

        // 4. Valores e vetores próprios
        double[][] valProATxA = valoresPropriosATxA(desviosA, desviosAT);
        double[][] vetProATxA = vetoresPropriosATxA(desviosA, desviosAT);
        print_Matrix(valProATxA, "Valores Próprios de A^T . A");
        print_Matrix(vetProATxA, "Vetores Próprios de A^T . A");

        double[][] vetProAxAT = vetoresPropriosAxAT(desviosA, desviosAT);
        double[][] vetNormalizados = normalizarVetores(vetProAxAT);
        print_Matrix(vetProAxAT, "Vetores Próprios de A . A^T");
        print_Matrix(vetNormalizados, "Vetores Próprios Normalizados");

        // 5. Centralização da imagem para reconstrução
        double[] linearizedOriginal = matrixToArray(arrayOriginal);
        System.out.println("Tamanho do linearizedOriginal: " + linearizedOriginal.length);
// Calcula o vetor médio para o vetor linearizado
        double[] meanVector = matrixToArray(arrayOriginal);
        System.out.println("Tamanho do vetor médio (meanVector): " + meanVector.length);


// Centraliza a imagem original linearizada
        double[] phi = centralizeImage(linearizedOriginal, meanVector);
        System.out.println("Tamanho de phi: " + phi.length);

        // 6. Cálculo dos pesos (projeção nos eigenfaces)
        double[] weights = calculateWeights(phi, vetNormalizados);
        System.out.println("Pesos Calculados (W): " + Arrays.toString(weights));

        // 7. Reconstrução da imagem com k eigenfaces
        double[] reconstructedImage = reconstructImage(meanVector, vetNormalizados, weights, k);
        System.out.println("Tamanho de reconstructedImage: " + reconstructedImage.length);
        print_Matrix(new double[][]{reconstructedImage}, "Imagem Reconstruída (1D)");

        // 8. Conversão para matriz 2D
        System.out.println("Tamanho de reconstructedImage: " + reconstructedImage.length);
        System.out.println("Tamanho esperado: " + (arrayOriginal.length * arrayOriginal[0].length));

        System.out.println("Dimensões da imagem original: " + arrayOriginal.length + "x" + arrayOriginal[0].length);
        System.out.println("Tamanho esperado: " + (arrayOriginal.length * arrayOriginal[0].length));
        System.out.println("Tamanho do vetor reconstruído: " + reconstructedImage.length);
        System.out.println("Tamanho do vetor médio (meanVector): " + meanVector.length);
        System.out.println("Dimensões de vetNormalizados: " + vetNormalizados.length + "x" + vetNormalizados[0].length);
        System.out.println("Tamanho dos pesos (weights): " + weights.length);


        double[][] reconstructedImageMatrix = arrayToMatrix(reconstructedImage, arrayOriginal.length, arrayOriginal[0].length);
        print_Matrix(reconstructedImageMatrix, "Imagem Reconstruída (2D)");
        print_Matrix(reconstructedImageMatrix, "Imagem Reconstruída (2D)");

        // 9. Salvar a imagem reconstruída
        saveImage(reconstructedImageMatrix, csvPath, "Input/TesteFuncao2-3/OutputImagesFunc2-3");
    }

    //* TESTE PARA SALVAR IMAGEM
    public static double[] matrixToArray(double[][] matrix) {
        int rows = matrix.length;
        int cols = matrix[0].length;
        double[] array = new double[rows * cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                array[i * cols + j] = matrix[i][j];
            }
        }
        return array;
    }


    public static double[] extractMeanVector(double[][] colunaMedia) {
        double[] meanVector = new double[colunaMedia[0].length];
        for (int i = 0; i < colunaMedia[0].length; i++) {
            meanVector[i] = colunaMedia[0][i]; // Extrai o primeiro elemento de cada coluna
        }
        return meanVector;
    }

    public static void saveImage(double[][] imageArray, String inputCsvPath, String outputFolderPath) {
        int height = imageArray.length;
        int width = imageArray[0].length;

        // Normalize the image to the range 0-255
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;

        // Find min and max values
        for (double[] row : imageArray) {
            for (double val : row) {
                if (val < min) min = val;
                if (val > max) max = val;
            }
        }

        int[][] normalizedImage = new int[height][width];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                normalizedImage[y][x] = (int) ((imageArray[y][x] - min) / (max - min) * 255);
            }
        }

        // Create output file name
        String csvFileName = new File(inputCsvPath).getName();
        String pngFileName = csvFileName.replace(".csv", ".png");
        String outputPath = outputFolderPath + "/" + pngFileName;

        // Write the image
        try {
            writeArrayAsImage(normalizedImage, outputPath);
            System.out.println("Image saved successfully: " + outputPath);
        } catch (IOException e) {
            System.err.println("Failed to save the image: " + e.getMessage());
        }
    }

    public static double[] centralizeImage(double[] image, double[] meanVector) {
        if (image.length != meanVector.length) {
            throw new IllegalArgumentException("Os tamanhos dos arrays 'image' e 'meanVector' não correspondem.");
        }
        double[] phi = new double[image.length];
        for (int i = 0; i < image.length; i++) {
            phi[i] = image[i] - meanVector[i];
        }
        return phi;
    }

    public static double[] calculateWeights(double[] phi, double[][] eigenfaces) {
        if (Math.sqrt(phi.length) != eigenfaces.length) {
            throw new IllegalArgumentException("O comprimento de 'phi' deve ser igual ao número de linhas em 'eigenfaces'.");
        }
        double[] weights = new double[eigenfaces[0].length];
        for (int j = 0; j < eigenfaces[0].length; j++) {
            for (int i = 0; i < eigenfaces[0].length; i++) {
                weights[j] += phi[i] * eigenfaces[i][j];
            }
        }
        return weights;
    }

    public static double[] reconstructImage(double[] meanVector, double[][] eigenfaces, double[] weights, int k) {
        // Limita k ao número de eigenfaces disponíveis
        k = Math.min(k, eigenfaces[0].length);

        // Inicializa o vetor reconstruído com o vetor médio
        double[] reconstructed = Arrays.copyOf(meanVector, meanVector.length);

        // Soma a contribuição de cada eigenface
        for (int j = 0; j < k; j++) {
            for (int i = 0; i < eigenfaces[0].length; i++) {
                reconstructed[i] += weights[j] * eigenfaces[i][j];
            }
        }

        return reconstructed;
    }




    public static double[][] arrayToMatrix(double[] array, int rows, int cols) {
        if (array.length != rows * cols) {
            throw new IllegalArgumentException("O tamanho do vetor não corresponde às dimensões da matriz.");
        }
        double[][] matrix = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                matrix[i][j] = array[i * cols + j];
            }
        }
        return matrix;
    }



    // Provided function to write an image
    public static void writeArrayAsImage(int[][] array, String outputFilePath) throws IOException {
        int height = array.length;
        int width = array[0].length;

        // Create a BufferedImage
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);

        // Set the pixel intensities
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int intensity = array[y][x];
                if (intensity < 0 || intensity > 255) {
                    throw new IllegalArgumentException("Pixel intensity must be between 0 and 255.");
                }
                int rgb = (intensity << 16) | (intensity << 8) | intensity; // Set the same value for R, G, B
                image.setRGB(x, y, rgb);
            }
        }

        // Write the image to the file
        File outputFile = new File(outputFilePath);
        ImageIO.write(image, "png", outputFile);
    }

    //* TESTE PARA SALVAR IMAGEM

    private static double[][] readCSVToArray(String filePath) {
        ArrayList<double[]> rows = new ArrayList<>();
        try (Scanner scanner = new Scanner(new File(filePath))) {
            while (scanner.hasNextLine()) {
                String line = scanner.nextLine().trim();
                if (!line.isEmpty()) {
                    String[] values = line.split(",");
                    double[] row = new double[values.length];
                    for (int i = 0; i < values.length; i++) {
                        row[i] = Double.parseDouble(values[i].trim());
                    }
                    rows.add(row);
                }
            }
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Erro: Arquivo CSV não encontrado no caminho: " + filePath, e);
        } catch (NumberFormatException e) {
            throw new RuntimeException("Erro: Valores inválidos encontrados no arquivo CSV.", e);
        }

        // Converte a lista de linhas em uma matriz
        double[][] matrix = new double[rows.size()][];
        for (int i = 0; i < rows.size(); i++) {
            matrix[i] = rows.get(i);
        }
        return matrix;
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