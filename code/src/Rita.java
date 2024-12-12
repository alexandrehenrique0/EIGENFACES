import java.io.File;
import java.io.FileNotFoundException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

//--------------2.1.2.2 Reconstrução/Representação de Imagens utilizando Eigenfaces-------------------------------------

public class Rita {
    public static void main(String[] args) throws IOException {
        String imageFolderLocation = "C:\\Users\\rita\\Documents\\lapr1-24-25_DAB_02\\Input\\TesteFuncao2-3\\csv";

        double[][] imageVectors = convertImagesToVectors(imageFolderLocation);

        for (int i = 0; i < imageVectors.length; i++) {
            System.out.printf("Imagem %d - Tamanho do vetor: %d", i + 1, imageVectors[i].length);
            System.out.println();
            System.out.print("Conteúdo do vetor: ");
            for (int j = 0; j < imageVectors[i].length; j++) {
                System.out.print(imageVectors[i][j] + " ");
            }

            System.out.println();
        }

        // Matriz M

        double[][] imageMatrix = buildImageMatrixM(imageVectors);

        System.out.println("\nMatriz de Imagens M (colunas são imagens):");
        for (int i = 0; i < imageMatrix.length; i++) {
            for (int j = 0; j < imageMatrix[i].length; j++) {
                System.out.print(imageMatrix[i][j] + " ");
            }
            System.out.println();
        }
    }

//--------1.Transformar cada imagem da base de imagens, com N imagens, para o formato vetor-----------------------------

    public static double[][] convertImagesToVectors(String imageFolderLocation) throws IOException {
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

        double[][] imageVectors = new double[csvFiles.length][];

        for (int i = 0; i < csvFiles.length; i++) {
            imageVectors[i] = loadImageFromCSV(csvFiles[i]);
        }

        return imageVectors;
    }


    public static double[] loadImageFromCSV(File imageFile) throws IOException {
        double[] imageVector = new double[0];

        if (!imageFile.exists() || !imageFile.isFile()) {
            throw new FileNotFoundException(String.format("Arquivo não encontrado: %s", imageFile.getName()));
        }

        BufferedReader reader = new BufferedReader(new FileReader(imageFile));
        String line;
        int rowCount = 0;
        int cols = 0;


        while ((line = reader.readLine()) != null) {
            if (rowCount == 0) {
                cols = line.split(",").length;
            }
            rowCount++;
        }

        imageVector = new double[cols * rowCount];


        reader.close();
        reader = new BufferedReader(new FileReader(imageFile));
        int index = 0;

        while ((line = reader.readLine()) != null) {
            String[] pixels = line.split(",");
            for (String pixel : pixels) {
                String trimmedPixel = pixel.trim();

                if (isValidDouble(trimmedPixel)) {
                    imageVector[index] = Double.parseDouble(trimmedPixel);
                } else {
                    System.out.printf("Valor inválido no arquivo %d na linha: %d", imageFile.getName(),line);
                    imageVector[index] = 0;
                }

                index++;
            }
        }
        reader.close();

        return imageVector;
    }

    public static boolean isValidDouble(String str) {
        try {
            Double.parseDouble(str);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }

    }

//------2.Construir uma matriz de imagens, M, em que cada coluna da matriz ´e uma imagem da base de imagens-------------

    public static double[][] buildImageMatrixM(double[][] imageVectors) {
        int numPixels = imageVectors[0].length;
        int numImages = imageVectors.length;

        double[][] imageMatrixM = new double[numPixels][numImages];

        for (int i = 0; i < numImages; i++) {
            for (int j = 0; j < numPixels; j++) {
                imageMatrixM[j][i] = imageVectors[i][j];
            }
        }

        return imageMatrixM;
    }
}