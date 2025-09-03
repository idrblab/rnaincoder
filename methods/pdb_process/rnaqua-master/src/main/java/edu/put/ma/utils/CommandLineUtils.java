package edu.put.ma.utils;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public final class CommandLineUtils {

    private static final int DEFAULT_LINE_WIDTH_OF_HELP = 200;
    private static final Logger LOGGER = LoggerFactory.getLogger(CommandLineUtils.class);

    private CommandLineUtils() {
        // hidden constructor
    }

    public static final void printHelp(final String artifactId, final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.setWidth(DEFAULT_LINE_WIDTH_OF_HELP);
        formatter.printHelp(artifactId, options);
        System.exit(1);
    }

    public static final CommandLine parseArgs(final String[] args, final Options options) {
        final CommandLineParser parser = new DefaultParser();
        CommandLine commandLine = null;
        try {
            commandLine = parser.parse(options, args);
        } catch (ParseException e) {
            LOGGER.error(e.getMessage(), e);
        }
        return commandLine;
    }
}
