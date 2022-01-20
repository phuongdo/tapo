package apps;

import org.apache.commons.cli.*;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.repeats.FinderMode;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Created by pdoviet on 10/10/2014.
 */
public class TaPo {
    static Logger logger = LoggerFactory.getLogger(TaPo.class);
    public static void main(String[] args) {
        String pdbCode = "1bpo";
        String pdbChain = "A";
        final CommandLineParser cmdLineGnuParser = new GnuParser();
        final Options options = constructOptions();
        CommandLine commandLine;
        String outputFile = null;
        int nCore = 0;
        try {
            commandLine = cmdLineGnuParser.parse(options, args);

            if (commandLine.hasOption("help")) {
                printUsage(options);
            } else {
                if (commandLine.hasOption("o")) {
                    outputFile = commandLine.getOptionValue("o");
                }
                if (commandLine.hasOption("nCore")) {
                    nCore = Integer.parseInt(commandLine
                            .getOptionValue("nCore"));

                }
                if (commandLine.hasOption("p")) {
                    pdbCode = commandLine.getOptionValue("p");
                    System.out.println(pdbCode);
                }
                if (commandLine.hasOption("c")) {
                    pdbChain = commandLine.getOptionValue("c");
                    System.out.println(commandLine.getOptionValue("c"));
                }

                if (commandLine.hasOption("f")) {
                    System.out.println(commandLine.getOptionValue("c"));
                }
                if (nCore == 0
                        || outputFile == null)
                    throw new org.apache.commons.cli.ParseException(
                            "Parse Exception!!!");
                RepeatFinder repeatFinder;
                if (commandLine.hasOption("f")) {
                    String fileDir = commandLine.getOptionValue("f");
//                    System.out.println("Loading file " + fileDir);
                    repeatFinder = new RepeatFinder(fileDir, pdbCode, pdbChain, nCore);
                    repeatFinder.setMode(FinderMode.CONSOLE);
                } else {
                    // load from PDF
//                    System.out.println("Loal PDB Bank");
                    repeatFinder = new RepeatFinder(pdbCode, pdbChain, nCore);
                    repeatFinder.setMode(FinderMode.CONSOLE);
                    repeatFinder.findRepeat();
                }
                String output;
                if (repeatFinder.isRepeat())
                    output = repeatFinder.getConsoleOutputDetails();
                else
                    output = pdbCode + "_" + pdbChain + "\tNo-TRs";
                DataIO.writeToFile(">" + pdbCode + "_" + pdbChain + "|" + repeatFinder.getTAPOScore() + "\n" + output, outputFile);
            }

        } catch (org.apache.commons.cli.ParseException e) {
            printUsage(options);
        } catch (Exception e) {
            // printUsage(options);
            e.printStackTrace();
        }

    }

    private static Options constructOptions() {
        final Options options = new Options();
        options.addOption("help", false, "help")
                .addOption("nCore", true, "number of cores you want to use. Default is 1 core")
                .addOption("p", true, "protein code eg. 1bpo")
                .addOption("c", true, "protein chain id eg. A")
                .addOption("f", true, "PDF file format")
                .addOption("o", true, "ouput file.");
        return options;
    }

    public static void printUsage(Options options) {
        HelpFormatter formatter = new HelpFormatter();
        String header = "available options as follow:\n\n";
        String footer = "\nPlease report issues to phuongdoviet@hotmail.com";
        formatter.printHelp("java -cp \"target/classes;lib/*\" apps.TaPo", header, options, footer,
                true);
        System.exit(1);
    }


}
