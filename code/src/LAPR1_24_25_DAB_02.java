import org.la4j.Matrix;
import java.io.File;
import java.util.Scanner;

public class LAPR1_24_25_DAB_02 {
    public static Scanner SCANNER = new Scanner(System.in);
    public static void main(String[] args) {
        int function = 0;
        int vectorNumbers = 0;
        String csvLocation = "";
        String imageLocation = "";

        if (check_Correct_Parameters(args)){
            System.out.println("Parâmetros corretos");
        }
        else if (args.length == 0)
        {
            // Mostrar as opções num menu e receber os parâmetros
            ui_Function_Parameter_Menu();
            function = receive_Function(null);
            
            ui_Vector_Numbers_Parameter_Menu();
            vectorNumbers = receive_Number_Vectors(null);
            
            ui_CSV_Location_Parameter_Menu();
            csvLocation = receive_CSV_Location(null);
            
            ui_Image_Location_Parameter_Menu();
            imageLocation = receive_Image_Location(null);
            // ---------------------------------------
            
        }else{
            System.out.println("Parâmetros inválidos");
            System.exit(1);
        }
    }

    //* Input --------------------------------------------------------------------------------
    //* --------------------------------------------------------------------------------------

    // Verificações de parâmetros
    public static boolean check_Correct_Parameters(String[] parameters){
        if (parameters.length == 8){
            return parameters[0].equals("-f") && parameters[2].equals("-k") && parameters[4].equals("-i") && parameters[6].equals("-j");
        }
        return false;
    }
    public static boolean check_function(int function){
        return function >= 1 && function <= 3;
    }
    public static boolean check_csvLocation(String csvLocation){
        File csv = new File(csvLocation);
        if (csvLocation.equals("")){
            return false;
        }else if (!csvLocation.contains(".csv")){
            return false;
        }
        else return csv.exists();
    }
    public static boolean check_imageLocation(String imageLocation){
        File imageDirectory = new File(imageLocation);
        if (imageLocation.isEmpty()){
            return false;
        }return imageDirectory.exists();
    }
    //--------------------------------------------------------

    // Menus de opções
    public static void  ui_Function_Parameter_Menu() {
        System.out.println("------------- Que função deseja realizar? -------------");
        System.out.println("1 - Decomposição Própria de uma Matriz Simétrica");
        System.out.println("2 - Reconstrução de Imagens usando Eigenfaces");
        System.out.println("3 - Identificação de imagem mais próxima");
        System.out.println("-------------------------------------------------------");
        System.out.printf("Opção: ");
    }
    public static void  ui_Vector_Numbers_Parameter_Menu() {
        System.out.println("------ Quantos vetores próprios deseja utilizar? ------");
        System.out.printf("Quantidade: ");
    }
    public static void  ui_CSV_Location_Parameter_Menu() {
        System.out.println("---- Qual a localização do csv que deseja utilizar? ---");
        System.out.printf("Localização: ");
    }
    public static void  ui_Image_Location_Parameter_Menu() {
        System.out.println("-------- Qual a localização da base de imagens? -------");
        System.out.printf("Localização: ");
    }
    //--------------------------------------------------------

    // Receber parâmetros
    public static int receive_Function(String[] args){
        int functionArgs = 0;
        if (args==null){
            int function = SCANNER.nextInt();
            if (!check_function(function)){
                System.out.println("Opção inválida");
                System.exit(1);
            }
            return function;
        }
        else{
            functionArgs = Integer.parseInt(args[1]);
            if (!check_function(functionArgs)){
                System.out.println("Opção inválida");
                System.exit(1);
            }
            return functionArgs;
        }

    }
    public static int receive_Number_Vectors(String[] args){
        int vectorNumbersArgs = 0;
        if (args==null){
            int vectorNumbers = SCANNER.nextInt();
            return vectorNumbers;
        }
        else{
            vectorNumbersArgs = Integer.parseInt(args[3]);
            return vectorNumbersArgs;
        }
    }
    public static String receive_CSV_Location(String[] args){
        String csvLocationArgs = "";
        if (args==null){
            String csvLocation = SCANNER.next();
            if (!check_csvLocation(csvLocation)){
                System.out.println("Localização inválida");
                System.exit(1);
            }
            return csvLocation;
        }
        else{
            csvLocationArgs = args[5];
            if (!check_csvLocation(csvLocationArgs)){
                System.out.println("Localização inválida");
                System.exit(1);
            }
            return csvLocationArgs;
        }
    }
    public static String receive_Image_Location(String[] args){
        String imageLocationArgs = "";
        if (args==null){
            String imageLocation = SCANNER.next();
            if (!check_imageLocation(imageLocation)){
                System.out.println("Localização inválida ou ficheiro não existe");
                System.exit(1);
            }
            return imageLocation;
        }
        else{
            imageLocationArgs = args[7];
            if (!check_imageLocation(imageLocationArgs)){
                System.out.println("Localização inválida ou ficheiro não existe");
                System.exit(1);
            }
            return imageLocationArgs;
        }
    }
    //--------------------------------------------------------

    //* --------------------------------------------------------------------------------------
    //* --------------------------------------------------------------------------------------
}
