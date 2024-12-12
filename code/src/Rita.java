import java.io.File;
import java.io.FileNotFoundException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

//2.1.2.2 Reconstrução/Representação de Imagens utilizando Eigenfaces

public class Rita {
    public static void main(String[] args) throws IOException {
        String imageFolderLocation = "C:\\Users\\rita\\Documents\\lapr1-24-25_DAB_02\\Input\\TesteFuncao2-3\\csv";  // Caminho para a pasta das imagens

        // Transformar imagens em vetores
        double[][] imageVectors = convertImagesToVectors(imageFolderLocation);

        // Exibir o tamanho do vetor para as primeiras imagens (apenas para verificar)
        for (int i = 0; i < imageVectors.length; i++) {
            System.out.printf("Imagem %d - Tamanho do vetor: %d", i + 1, imageVectors[i].length);
            System.out.println();
            System.out.print("Conteúdo do vetor: ");
            for (int j = 0; j < imageVectors[i].length; j++) {
                System.out.print(imageVectors[i][j] + " ");
            }

            System.out.println(); // Quebra de linha para separar as imagens
        }

        // Construir a matriz de imagens M
        double[][] imageMatrix = buildImageMatrixM(imageVectors);

        // Exibir a matriz de imagens (só para verificar)
        System.out.println("\nMatriz de Imagens M (colunas são imagens):");
        for (int i = 0; i < imageMatrix.length; i++) {
            for (int j = 0; j < imageMatrix[i].length; j++) {
                System.out.print(imageMatrix[i][j] + " ");
            }
            System.out.println();
        }
    }

//--------1.Transformar cada imagem da base de imagens, com N imagens, para o formato vetor-----------------------------

    // Método para converter as imagens em vetores, ajustando o tamanho das imagens
    public static double[][] convertImagesToVectors(String imageFolderLocation) throws IOException {
        File folder = new File(imageFolderLocation);

        // Obter todos os arquivos no diretório
        File[] imageFiles = folder.listFiles();

        if (imageFiles == null || imageFiles.length == 0) {
            throw new FileNotFoundException("Nenhuma imagem CSV encontrada na pasta.");
        }

        // Contar quantos arquivos CSV existem
        int csvFileCount = 0;
        for (File file : imageFiles) {
            if (file.isFile() && file.getName().toLowerCase().endsWith(".csv")) {
                csvFileCount++;
            }
        }

        // Se nenhum arquivo CSV foi encontrado
        if (csvFileCount == 0) {
            throw new FileNotFoundException("Nenhum arquivo CSV encontrado.");
        }

        // Inicializar o array de arquivos CSV com o tamanho correto
        File[] csvFiles = new File[csvFileCount];
        int index = 0;

        // Adicionar os arquivos CSV ao array
        for (File file : imageFiles) {
            if (file.isFile() && file.getName().toLowerCase().endsWith(".csv")) {
                csvFiles[index++] = file;
            }
        }

        // Exibir quantidade de arquivos CSV encontrados
        System.out.printf("Encontrados %d arquivos CSV.%n", csvFiles.length);

        // Inicializar o array para armazenar todos os vetores das imagens
        double[][] imageVectors = new double[csvFiles.length][];

        // Para cada arquivo CSV, carregamos os dados e transformamos em vetores
        for (int i = 0; i < csvFiles.length; i++) {
            imageVectors[i] = loadImageFromCSV(csvFiles[i]);
        }

        return imageVectors;
    }


    // Carregar uma imagem de um arquivo CSV e transformá-la em vetor
    public static double[] loadImageFromCSV(File imageFile) throws IOException {
        double[] imageVector = new double[0];  // Inicializa como vetor vazio

        // Verificar se o arquivo existe antes de tentar ler
        if (!imageFile.exists() || !imageFile.isFile()) {
            throw new FileNotFoundException(String.format("Arquivo não encontrado: %s", imageFile.getName()));
        }


        // Lendo o arquivo CSV
        BufferedReader reader = new BufferedReader(new FileReader(imageFile));
        String line;
        int rowCount = 0;
        int cols = 0;

        // Primeiro, contar o número de colunas e linhas
        while ((line = reader.readLine()) != null) {
            if (rowCount == 0) {
                cols = line.split(",").length;  // número de colunas
            }
            rowCount++;
        }

        // Inicializar o vetor com o número correto de elementos (linhas * colunas)
        imageVector = new double[cols * rowCount];

        // Reposicionar o leitor para a leitura real dos dados
        reader.close();
        reader = new BufferedReader(new FileReader(imageFile));
        int index = 0;

        // Lendo os dados reais
        while ((line = reader.readLine()) != null) {
            String[] pixels = line.split(",");
            for (String pixel : pixels) {
                String trimmedPixel = pixel.trim();  // Remover espaços extras

                // Verificar se a conversão é possível
                if (isValidDouble(trimmedPixel)) {
                    imageVector[index] = Double.parseDouble(trimmedPixel);  // Atribuir o valor convertido
                } else {
                    System.out.printf("Valor inválido no arquivo %d na linha: %d", imageFile.getName(),line);
                    imageVector[index] = 0;  // Atribuir um valor padrão em caso de erro
                }

                index++;
            }
        }
        reader.close();

        return imageVector;
    }

    // Método para verificar se uma string pode ser convertida para double
    public static boolean isValidDouble(String str) {
        try {
            Double.parseDouble(str);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }

    }

//------2.Construir uma matriz de imagens, M, em que cada coluna da matriz ´e uma imagem da base de imagens-------------

    // Método para construir a matriz de imagens, onde cada coluna é uma imagem (vetor)
    public static double[][] buildImageMatrixM(double[][] imageVectors) {
        // Determinar o número de pixels de cada imagem e o número de imagens
        int numPixels = imageVectors[0].length;  // Número de pixels por imagem (tamanho do vetor de cada imagem)
        int numImages = imageVectors.length;     // Número de imagens

        // Inicializar a matriz de imagens M com tamanho numPixels x numImages
        double[][] imageMatrixM = new double[numPixels][numImages];

        // Preencher a matriz com os vetores das imagens
        for (int i = 0; i < numImages; i++) {
            for (int j = 0; j < numPixels; j++) {
                imageMatrixM[j][i] = imageVectors[i][j];  // A coluna i é a imagem i
            }
        }

        return imageMatrixM;
    }
}