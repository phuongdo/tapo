package apps;

import org.apache.commons.cli.*;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.parallel.ExecutorServiceApp;
import org.cnrs.crbm.lib.parallel.JobInput;
import org.cnrs.crbm.lib.parallel.JobOuput;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;

public class AnalysisApp {

	static Logger logger = LoggerFactory.getLogger(AnalysisApp.class);
	public static void main(String[] args) {

		// Set up a simple configuration that logs on the console.
		// PropertiesConfigurator is used to configure logger from properties
		// file
		//
		// get a logger instance named "com.foo"
		// Now set its level. Normally you do not need to set the
		// level of a logger programmatically. This is usually done
		// in configuration files.

		final CommandLineParser cmdLineGnuParser = new GnuParser();
		final Options options = constructOptions();
		CommandLine commandLine;
		String inputFile = null, outputFile = null;
		int nJob = 0;
		int nCore = 0;

		try {
			commandLine = cmdLineGnuParser.parse(options, args);

			if (commandLine.hasOption("help")) {
				printUsage(options);
			} else {

				if (commandLine.hasOption("i")) {
					inputFile = commandLine.getOptionValue("i");

				}
				if (commandLine.hasOption("o")) {
					outputFile = commandLine.getOptionValue("o");

				}
				if (commandLine.hasOption("nJob")) {
					nJob = Integer.parseInt(commandLine.getOptionValue("nJob"));

				}
				if (commandLine.hasOption("nCore")) {
					nCore = Integer.parseInt(commandLine
							.getOptionValue("nCore"));

				}

				if (nJob == 0 || nCore == 0 || inputFile == null
						|| outputFile == null)
					throw new org.apache.commons.cli.ParseException(
							"Parse Exception!!!");

				logger.info("Starting program...");
				logger.info("parameters: nJob:" + nJob + " nCore:" + nCore);
				// your program here

				logger.info("Download your data...");

				// ProteinCSVReader csvReader = new ProteinCSVReader();

				// List<Row> rows = csvReader.getData(inputFile);

				// ReadFasta fastaReader = new ReadFasta();
				// List<Row> rows = fastaReader
				// .getRowAllPdbsWithChain(Dir.FASTA_LOCAL
				// + "/pdb_seqres.nrdb");

				List<Row> rows = DataIO.getListProteinFromFile(inputFile);
				JobInput jobInput = new JobInput();
				jobInput.setProteins(rows);

				logger.info("Processing...");
				ExecutorServiceApp executor = new ExecutorServiceApp(nJob,
						nCore, jobInput);
				JobOuput joboutput = executor.run();

				// logger.info(joinedOutput);
				StringBuilder builder = new StringBuilder();
				// builder.append("<table>");
				// builder.append("<tr><td>protein</td><td>Tag</td><td>Conclusion</td><td>Annotation</td><td>TRs</td><td>RMSD</td><td>T-REKs 3D</td><td>TRUST 3D</td><td>TRUST 3D(*)</td><td>VECTORS</td><td>VECTORS-2ndS</td><td>CONTACT MAP</td><td>RAPHEAL</td><tr/>");
				// builder.append("protein\tTRs\tanno\tlength\n");
				builder.append(joboutput);
				// builder.append("</table>");

				// DataIO.writeToFile(builder.toString(), outputFile);
				// logger.info("Stop program....");
				// System.exit(1);
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
				.addOption("nCore", true, "number of cores you want to use.")
				.addOption("nJob", true,
						"number of jobs that were splitted from list of proteins")
				.addOption("o", true, "output of program")
				.addOption("i", true, "input list of pdbs files.");
		return options;
	}

	public static void printUsage(Options options) {
		HelpFormatter formatter = new HelpFormatter();

		String header = "available options as follow:\n\n";
		String footer = "\nPlease report issues to phuongdoviet@hotmail.com";
		formatter.printHelp("java -jar StructLib.jar", header, options, footer,
				true);
		System.exit(1);
	}
}
