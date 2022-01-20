package apps;

import org.apache.commons.cli.*;
import org.junit.Test;

public class GnuParserTest {





    @Test
    public void testParser() throws ParseException {
        final CommandLineParser cmdLineGnuParser = new GnuParser();
        final Options options = constructOptions();
        CommandLine commandLine;
        String[] args="-p 1ib2 -c A -nCore 1 -o /home/sgeadmin/work/1ib2.o".split(" ");
        commandLine = cmdLineGnuParser.parse(options, args);
        System.out.println(commandLine.getOptionValue("o"));

    }

    private static Options constructOptions() {
        final Options options = new Options();
        options.addOption("help", false, "help")
                .addOption("nCore", true, "number of cores you want to use. Default is 1 core")
                .addOption("p", true, "protein code eg. 1bpo")
                .addOption("c", true, "protein chain id eg. A")
//                .addOption("f", false, "PDF file format")
                .addOption("o", true, "ouput file.");
        return options;
    }
}
