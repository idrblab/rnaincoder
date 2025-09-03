package edu.put.ma.model.input;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.collect.ImmutableSet;

import edu.put.ma.Command;
import edu.put.ma.CommandEnum;
import edu.put.ma.utils.ArgumentUtils;
import edu.put.ma.utils.CommandLineUtils;
import edu.put.ma.utils.PreconditionUtils;

public final class InputModelsFactory {

    private static final Logger LOGGER = LoggerFactory.getLogger(InputModelsFactory.class);

    private static final ImmutableSet<String> COMMAND_CODES = ImmutableSet.of("c", "command");

    private static final ImmutableSet<String> REFERENCE_STRUCTURE_CODES = ImmutableSet.of("r",
            "reference-structure-file-path");

    private static final ImmutableSet<String> SINGLE_MODEL_CODES = ImmutableSet.of("s",
            "single-model-file-path");

    private static final ImmutableSet<String> MULTIPLE_MODELS_CODES = ImmutableSet.of("d",
            "multiple-models-directory-path");

    private InputModelsFactory() {
        // hidden constructor
    }

    public static final CommonInputModel getInputModel(final String[] args, final String artifactId) {
        final Command command = getCommand(args, artifactId);
        if (command == CommandEnum.DEFORMATION_PROFILE) {
            return new DeformationProfileInputModelImpl(args, artifactId);
        } else {
            final String[] singleModelArg = ArgumentUtils.retrieveArgByOpt(args, SINGLE_MODEL_CODES);
            final String[] multipleModelsArg = ArgumentUtils.retrieveArgByOpt(args, MULTIPLE_MODELS_CODES);
            final String[] referenceStructureArg = ArgumentUtils.retrieveArgByOpt(args,
                    REFERENCE_STRUCTURE_CODES);
            if (!CollectionUtils.sizeIsEmpty(referenceStructureArg)) {
                return new ComparisonInputModelImpl(args, artifactId);
            } else if ((!CollectionUtils.sizeIsEmpty(singleModelArg))
                    || (!CollectionUtils.sizeIsEmpty(multipleModelsArg))) {
                return new AnalysisInputModelImpl(args, artifactId);
            }
        }
        throw new IllegalArgumentException("Inappropriate command definition!");
    }

    private static final Options getCommandOption() {
        final Options options = new Options();
        options.addOption("c", "command", true, CommonInputModelImpl.getCommandDescriptionString());
        return options;
    }

    private static final Command getCommand(final Options options, final String[] commandArg,
            final String artifactId) {
        final String commandCode = getCommandCode(commandArg[0]);
        final CommandLine commandLine = CommandLineUtils.parseArgs(commandArg, options);
        if (commandLine.hasOption(commandCode)) {
            final Command command = CommonInputModelImpl.getCommand(commandLine, commandCode);
            if (command != null) {
                return command;
            } else {
                CommandLineUtils.printHelp(artifactId, options);
            }
        }
        throw new IllegalArgumentException("No command defined!");
    }

    private static final String getCommandCode(final String code) {
        PreconditionUtils.checkIfStringIsBlank(code, "Command");
        final String longPrefix = "--";
        final String shortPrefix = "-";
        String result = !StringUtils.startsWith(code, longPrefix) ? code : StringUtils.substring(code,
                StringUtils.length(longPrefix));
        result = !StringUtils.startsWith(result, shortPrefix) ? result : StringUtils.substring(result,
                StringUtils.length(shortPrefix));
        return result;
    }

    private static final Command getCommand(final String[] args, final String artifactId) {
        final Options options = getCommandOption();
        try {
            final String[] commandArg = ArgumentUtils.retrieveArgByOpt(args, COMMAND_CODES);
            if (!CollectionUtils.sizeIsEmpty(commandArg)) {
                return getCommand(options, commandArg, artifactId);
            } else {
                CommandLineUtils.printHelp(artifactId, options);
            }
        } catch (Exception e) {
            LOGGER.error(e.getMessage(), e);
            CommandLineUtils.printHelp(artifactId, options);
        }
        throw new IllegalArgumentException("No command defined!");
    }

}
