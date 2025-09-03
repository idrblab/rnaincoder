package edu.put.ma.model.input;

import static edu.put.ma.utils.StringUtils.NEW_LINE;

import java.io.File;
import java.lang.reflect.Field;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import edu.put.ma.Command;
import edu.put.ma.CommandEnum;
import edu.put.ma.model.services.Fragment;
import edu.put.ma.model.services.PdbStructure;
import edu.put.ma.model.services.StructuresSet;
import edu.put.ma.utils.ArrayUtils;
import edu.put.ma.utils.CommandLineUtils;
import edu.put.ma.utils.PreconditionUtils;

public abstract class CommonInputModelImpl implements CommonInputModel {

    private static final Logger LOGGER = LoggerFactory.getLogger(CommonInputModelImpl.class);

    private Options options;

    private String targetModelsAlignment;

    private CommandEnum.BpsIdentificationTool basePairsIdentificationTool;

    private boolean considerAtomsSupportedByRNAPuzzlesOnly;

    private CommandLine commandLine;

    @Getter
    private Command command;

    @Getter
    @Setter
    private String outputFilePath;

    protected static final int PAIR_SIZE = 2;

    protected Map<String, String> optionsMapping;

    protected CommonInputModelImpl(final String[] args) {
        this.options = constructCommonOptions();
        extendOptions(options, constructSpecificOptions());
        this.commandLine = CommandLineUtils.parseArgs(args, options);
    }

    protected CommonInputModelImpl(final Builder commonInputModelBuilder) {
        this.command = commonInputModelBuilder.command;
        this.outputFilePath = commonInputModelBuilder.outputFilePath;
        this.targetModelsAlignment = commonInputModelBuilder.targetModelsAlignment;
        this.basePairsIdentificationTool = commonInputModelBuilder.basePairsIdentificationTool;
        this.considerAtomsSupportedByRNAPuzzlesOnly = commonInputModelBuilder.considerAtomsSupportedByRNAPuzzlesOnly;
        this.initOptionsMapping();
    }

    @Override
    public StructuresSet getStructuresSet() {
        final StructuresSet structuresSet = constructStructuresSet();
        introduceTargetVsModelsAlignment(structuresSet);
        return structuresSet;
    }

    @Override
    public StructuresSet constructStructuresSet() {
        final StructuresSet structuresSet = new StructuresSet();
        structuresSet.setAnnotationTool(basePairsIdentificationTool.getMapping().getServiceCode());
        structuresSet.setIsFiltering(considerAtomsSupportedByRNAPuzzlesOnly);
        return structuresSet;
    }

    @Override
    public void printHelp(final String artifactId) {
        CommandLineUtils.printHelp(artifactId, options);
    }

    @Override
    public String[] getArgs() {
        final List<String> result = Lists.newArrayList();
        Class<?> currentClass = this.getClass();
        while (currentClass.getSuperclass() != null) {
            CollectionUtils.addAll(result, processDeclaredFields(currentClass, this));
            currentClass = currentClass.getSuperclass();
        }
        return result.toArray(new String[CollectionUtils.size(result)]);
    }

    @Override
    public boolean isInputInitializedProperly() {
        return (isCommandLineHasOption("c")) && (isCommandLineHasOption("o"))
                && isTargetVsModelsAlignmentValid() && isOutputPathValid();
    }

    @Override
    public String getServiceCode() {
        return command.getMapping().getServiceCode();
    }

    @Override
    public String getCommandName() {
        return command.getMapping().getName();
    }

    @Override
    public String getInputModelString() {
        final StringBuilder inputModelStringBuilder = new StringBuilder("Command: ").append(command);
        edu.put.ma.utils.StringUtils.extendStringBuilder(inputModelStringBuilder, targetModelsAlignment,
                "Alignment");
        if (basePairsIdentificationTool != null) {
            inputModelStringBuilder.append("\nBase pairs identification tool: ").append(
                    basePairsIdentificationTool);
        }
        inputModelStringBuilder.append("\nConsider atoms supported by RNA-Puzzles only: ").append(
                considerAtomsSupportedByRNAPuzzlesOnly ? "Y" : "N");
        edu.put.ma.utils.StringUtils.extendStringBuilder(inputModelStringBuilder, outputFilePath,
                "Output file path");
        return inputModelStringBuilder.append(NEW_LINE).toString();
    }

    protected void introduceTargetVsModelsAlignment(final StructuresSet structuresSet) {
        if (structuresSet.getStructureListSize() > 0) {
            final Map<String, List<Fragment>> targetVsModelsAlignment = this.getTargetVsModelsAlignment();
            if (!CollectionUtils.sizeIsEmpty(targetVsModelsAlignment)) {
                for (PdbStructure structure : structuresSet.getStructuresList()) {
                    setFragmentList(targetVsModelsAlignment, structure);
                }
            }
        }
    }

    private void setFragmentList(final Map<String, List<Fragment>> targetVsModelsAlignment,
            final PdbStructure structure) {
        final String structureFilename = structure.getFilename();
        if (targetVsModelsAlignment.containsKey(structureFilename)) {
            structure.setFragmentList(targetVsModelsAlignment.get(structureFilename));
        }
    }

    protected boolean isOptionSet(final String option) {
        return isOptionSet(commandLine, option);
    }

    protected Map<String, List<Fragment>> getTargetVsModelsAlignment() {
        if (StringUtils.isNotBlank(targetModelsAlignment)) {
            return getTargetVsModelsAlignment(targetModelsAlignment);
        }
        return Collections.emptyMap();
    }

    protected void initOptionsMapping() {
        this.optionsMapping = Maps.newHashMap();
        this.optionsMapping.putAll(ImmutableMap.of("command", "-c", "outputFilePath", "-o",
                "targetModelsAlignment", "-a", "basePairsIdentificationTool", "-t",
                "considerAtomsSupportedByRNAPuzzlesOnly", "-f"));
    }

    protected boolean isCommandLineHasOption(final String option) {
        return commandLine.hasOption(option);
    }

    protected String getOptionString(final String option) {
        return getOptionString(commandLine, option);
    }

    protected void secureInitState(final String artifactId) {
        try {
            initState();
        } catch (Exception e) {
            LOGGER.error(e.getMessage(), e);
            printHelp(artifactId);
        }
    }

    protected void initState() {
        setCommand();
        setBasePairsIdentificationTools();
        setTargetModelsAlignment();
        setConsiderAtomsSupportedByRNAPuzzlesOnly();
        setOutputFilePath();
    }

    @Getter
    @NoArgsConstructor
    public static class Builder {

        private Command command;

        private String outputFilePath;

        private String targetModelsAlignment;

        private CommandEnum.BpsIdentificationTool basePairsIdentificationTool;

        private boolean considerAtomsSupportedByRNAPuzzlesOnly;

        public Builder command(final Command command) {
            this.command = command;
            return this;
        }

        public Builder outputFilePath(final String outputFilePath) {
            this.outputFilePath = outputFilePath;
            return this;
        }

        public Builder targetModelsAlignment(final String targetModelsAlignment) {
            this.targetModelsAlignment = targetModelsAlignment;
            return this;
        }

        public Builder basePairsIdentificationTool(
                final CommandEnum.BpsIdentificationTool basePairsIdentificationTool) {
            this.basePairsIdentificationTool = basePairsIdentificationTool;
            return this;
        }

        public Builder considerAtomsSupportedByRNAPuzzlesOnly(
                final boolean considerAtomsSupportedByRNAPuzzlesOnly) {
            this.considerAtomsSupportedByRNAPuzzlesOnly = considerAtomsSupportedByRNAPuzzlesOnly;
            return this;
        }
    }

    private List<String> processDeclaredFields(final Class<?> currentClass, final Object object) {
        final List<String> result = Lists.newArrayList();
        for (Field field : currentClass.getDeclaredFields()) {
            field.setAccessible(true);
            try {
                final String fieldName = field.getName();
                final Object fieldObject = field.get(object);
                if ((this.optionsMapping.containsKey(fieldName))
                        && (isConsideredString(field, fieldObject)
                                || isConsideredPrimitive(field, fieldObject) || isNotNullOtherInstance(field,
                                    fieldObject))) {
                    addParameter(result, field, fieldName, fieldObject);
                }
            } catch (IllegalArgumentException e) {
                LOGGER.error(e.getMessage(), e);
            } catch (IllegalAccessException e) {
                LOGGER.error(e.getMessage(), e);
            }
        }
        return result;
    }

    private void addParameter(final List<String> result, final Field field, final String fieldName,
            final Object fieldObject) {
        if (!isConsideredBoolean(field)) {
            result.add(optionsMapping.get(fieldName));
            result.add(String.valueOf(fieldObject));
        } else if (Boolean.valueOf(String.valueOf(fieldObject)).booleanValue()) {
            result.add(optionsMapping.get(fieldName));
        }
    }

    private static boolean isNotNullOtherInstance(final Field field, final Object fieldObject) {
        return (!((field.getType().isAssignableFrom((Class<?>) Integer.TYPE)) || (field.getType()
                .isAssignableFrom((Class<?>) Double.TYPE)))) && (fieldObject != null);
    }

    private static boolean isConsideredString(final Field field, final Object fieldObject) {
        return (field.getType().isAssignableFrom((Class<?>) String.class)) && (fieldObject != null)
                && (StringUtils.isNotBlank(String.valueOf(fieldObject)));
    }

    private static boolean isConsideredPrimitive(final Field field, final Object fieldObject) {
        return ((field.getType().isAssignableFrom((Class<?>) Integer.TYPE)) && (((Integer) fieldObject)
                .intValue() > 0))
                || ((field.getType().isAssignableFrom((Class<?>) Double.TYPE)) && (Double.compare(
                        ((Double) fieldObject).doubleValue(), 0.0) > 0));
    }

    private static boolean isConsideredBoolean(final Field field) {
        return field.getType().isAssignableFrom((Class<?>) Boolean.TYPE);
    }

    private void setCommand() {
        command = getCommand(commandLine, "c");
    }

    private void setTargetModelsAlignment() {
        targetModelsAlignment = getOptionString("a");
    }

    private void setBasePairsIdentificationTools() {
        basePairsIdentificationTool = getBasePairsIdentificationTool(commandLine, "t");
    }

    private void setConsiderAtomsSupportedByRNAPuzzlesOnly() {
        considerAtomsSupportedByRNAPuzzlesOnly = isOptionSet(commandLine, "f");
    }

    private void setOutputFilePath() {
        outputFilePath = getOptionString("o");
    }

    private void extendOptions(Options options, Options additionalOptions) {
        for (Option option : additionalOptions.getOptions()) {
            options.addOption(option);
        }
    }

    private boolean isOutputPathValid() {
        final File outputFile = FileUtils.getFile(outputFilePath);
        PreconditionUtils.checkIfFileDoesNotExistOrIsNotADirectory(outputFile, "Output file");
        return true;
    }

    private boolean isTargetVsModelsAlignmentValid() {
        if (StringUtils.isNotBlank(targetModelsAlignment)) {
            final String currentAlign = edu.put.ma.utils.StringUtils.deleteDanglingChars(
                    targetModelsAlignment, ";");
            final String[] models = StringUtils.split(currentAlign, ';');
            boolean isTarget = true;
            final List<Integer> distances = Lists.newArrayList();
            final StringBuilder sb = new StringBuilder();
            for (String modelAlignment : models) {
                if (!isTarget) {
                    sb.append(';');
                }
                modelAlignment = edu.put.ma.utils.StringUtils.deleteDanglingChars(modelAlignment, "|");
                analyseModelAlignment(isTarget, distances, sb, modelAlignment);
                isTarget = false;
            }
            final String validatedAlignmentString = sb.toString();
            if (StringUtils.isNotBlank(validatedAlignmentString)) {
                this.targetModelsAlignment = validatedAlignmentString;
            } else {
                return false;
            }
        }
        return true;
    }

    public static final String getCommandDescriptionString() {
        return new StringBuilder("supported commands: ")
                .append(ArrayUtils.getEnumNamesString(CommandEnum.class)).append(", ")
                .append(ArrayUtils.getEnumNamesString(CommandEnum.Measure.class)).append(", ")
                .append(ArrayUtils.getEnumNamesString(CommandEnum.Sequence.class)).append(", ")
                .append(ArrayUtils.getEnumNamesString(CommandEnum.Structure3d.class)).toString();
    }

    public static final Command getCommand(final CommandLine commandLine, final String option) {
        Command command = getEnumName(commandLine, option, CommandEnum.class);
        command = (command == null) ? getEnumName(commandLine, option, CommandEnum.Measure.class) : command;
        command = (command == null) ? getEnumName(commandLine, option, CommandEnum.Sequence.class) : command;
        command = (command == null) ? getEnumName(commandLine, option, CommandEnum.Structure3d.class)
                : command;
        if (command != null) {
            return command;
        }
        throw new IllegalArgumentException(String.format("No command %s found!",
                getOptionString(commandLine, option)));
    }

    public static final boolean isOptionSet(final CommandLine commandLine, final String option) {
        return commandLine.hasOption(option);
    }

    private static <T extends Enum<T>> T getEnumValue(final CommandLine commandLine, final String option,
            final Class<T> enumClass) {
        return getEnumNameOrValue(commandLine, option, enumClass, false);
    }

    private static <T extends Enum<T>> T getEnumName(final CommandLine commandLine, final String option,
            final Class<T> enumClass) {
        return getEnumNameOrValue(commandLine, option, enumClass, true);
    }

    private static <T extends Enum<T>> T getEnumNameOrValue(final CommandLine commandLine,
            final String option, final Class<T> enumClass, final boolean isName) {
        try {
            final String optionString = getOptionString(commandLine, option);
            if (StringUtils.isNotBlank(optionString)) {
                T enumR = (isName) ? CommandEnum.fromString(optionString, enumClass) : getEnumValue(
                        optionString, enumClass);
                if (enumR != null) {
                    return enumR;
                }
            }
        } catch (Exception e) {
            LOGGER.debug(e.getMessage(), e);
        }
        return null;
    }

    private static void analyseModelAlignment(final boolean isTarget, final List<Integer> distances,
            final StringBuilder sb, final String modelAlignment) {
        final String[] modelAlignmentEntries = StringUtils.split(modelAlignment, ':');
        if (org.apache.commons.lang3.ArrayUtils.getLength(modelAlignmentEntries) == PAIR_SIZE) {
            final String modelName = StringUtils.trim(modelAlignmentEntries[0]);
            if (StringUtils.endsWithIgnoreCase(modelName, ".pdb")) {
                sb.append(modelName).append(':');
            } else {
                sb.append(modelName).append(".pdb:");
            }
            final String[] fragments = StringUtils.split(StringUtils.trim(modelAlignmentEntries[1]), '|');
            analyseAlignmentFragments(isTarget, distances, sb, fragments);
        } else {
            LOGGER.info("You should define the model alignment according to the following example 'model.pdb:A_1,5|A_10,5'.");
        }
    }

    private static void analyseAlignmentFragments(final boolean isTarget, final List<Integer> distances,
            final StringBuilder sb, final String[] fragments) {
        int fragmentIdx = 0;
        for (String fragment : fragments) {
            final String[] fragmentEntries = StringUtils.split(fragment, ',');
            final String distance = StringUtils.trim(fragmentEntries[1]);
            if (!Pattern.matches("[0-9]+", distance)) {
                LOGGER.info(String.format("Nucleotides distance (%d) should be integer number e.g., 5.",
                        distance));
            } else {
                parseFragment(isTarget, distances, sb, fragmentIdx, fragmentEntries, distance);
                fragmentIdx++;
            }

        }
    }

    private static void parseFragment(final boolean isTarget, final List<Integer> distances,
            final StringBuilder sb, final int fragmentIdx, final String[] fragmentEntries,
            final String distance) {
        boolean isOk = true;
        final String[] residueKeyEntries = StringUtils.split(fragmentEntries[0], '_');
        String chain = StringUtils.trim(residueKeyEntries[0]).replaceAll(" +", " ");
        chain = "".equals(chain) ? " " : chain;
        final String resNo = StringUtils.trim(residueKeyEntries[1]);
        final String iCode = (org.apache.commons.lang3.ArrayUtils.getLength(residueKeyEntries) >= 3) ? residueKeyEntries[2]
                .replaceAll(" +", " ") : " ";
        if (org.apache.commons.lang3.ArrayUtils.getLength(residueKeyEntries) >= PAIR_SIZE) {
            isOk = verifyChain(chain) && verifyResNo(resNo) && verifyIcode(iCode);
        } else {
            LOGGER.info("The first residue of each considered fragment should be described by chain + '_' + residue number (+ '_' + insertion code optionally).");
        }
        if (isOk) {
            if (isTarget) {
                distances.add(Integer.parseInt(distance));
                addFragment(sb, fragmentIdx, distance, chain, resNo, iCode);
            } else {
                if (Integer.parseInt(distance) == distances.get(fragmentIdx)) {
                    addFragment(sb, fragmentIdx, distance, chain, resNo, iCode);
                } else {
                    LOGGER.info(
                            String.format("Nucleotides distance defined for model (%d) does not correspond to the particular target distance (%d)."),
                            distances, distances.get(fragmentIdx));
                }
            }
        }
    }

    private static void addFragment(final StringBuilder sb, final int fragmentIdx, final String distance,
            final String chain, final String resNo, final String iCode) {
        if (fragmentIdx > 0) {
            sb.append('|');
        }
        sb.append(chain).append('_').append(resNo).append('_').append(iCode).append(',').append(distance);
    }

    private static boolean verifyIcode(final String iCode) {
        if (!Pattern.matches("[0-9A-Za-z ]", iCode)) {
            LOGGER.info(String.format(
                    "iCode (%s) should be a single-letter e.g., digit, space, lower- or upper-case letter.",
                    iCode));
            return false;
        }
        return true;
    }

    private static boolean verifyResNo(final String resNo) {
        if (!Pattern.matches("[-]?[0-9]+", resNo)) {
            LOGGER.info(String.format("ResNo (%d) should be an integer number e.g., 5.", resNo));
            return false;
        }
        return true;
    }

    private static boolean verifyChain(final String chain) {
        if (!Pattern.matches("[0-9A-Za-z ]", chain)) {
            LOGGER.info(String
                    .format("Chain (%s) should be a single-letter code e.g., digit, space, lower- or upper-case letter.",
                            chain));
            return false;
        }
        return true;
    }

    private static String getOptionString(final CommandLine commandLine, final String option) {
        if (commandLine.hasOption(option)) {
            return commandLine.getOptionValue(option);
        }
        return null;
    }

    private static Options constructCommonOptions() {
        final Options options = new Options();
        options.addOption("c", "command", true, getCommandDescriptionString());
        options.addOption("a", "alignment", true, "(optional) single- or multi-models alignment");
        options.addOption(
                "t",
                "base-pairs-identification-tool",
                true,
                new StringBuilder("(optional) base pairs identification tool, supported tools: ").append(
                        ArrayUtils.getEnumNamesString(CommandEnum.BpsIdentificationTool.class)).toString());
        options.addOption("f", "consider-atoms-supported-by-RNA-Puzzles-only", false,
                "(optional) consider atoms supported by RNA-Puzzles only");
        options.addOption("o", "output-file-path", true, "output file path");
        return options;
    }

    private static CommandEnum.BpsIdentificationTool getBasePairsIdentificationTool(
            final CommandLine commandLine, final String option) {
        final CommandEnum.BpsIdentificationTool bpsIdentificationTool = getEnumValue(commandLine, option,
                CommandEnum.BpsIdentificationTool.class);
        if (bpsIdentificationTool != null) {
            return bpsIdentificationTool;
        }
        return CommandEnum.BpsIdentificationTool.MC_ANNOTATE;
    }

    private static <T extends Enum<T>> T getEnumValue(final String value, final Class<T> enumClass) {
        if ((enumClass != null) && (StringUtils.isNotBlank(value))) {
            return Enum.valueOf(enumClass, StringUtils.upperCase(StringUtils.trim(value)));
        }
        return null;
    }

    private static Map<String, List<Fragment>> getTargetVsModelsAlignment(final String targetVsModelsAlignment) {
        final Map<String, List<Fragment>> result = Maps.newHashMap();
        final String[] modelAlignments = StringUtils.split(targetVsModelsAlignment, ';');
        for (String modelAlignment : modelAlignments) {
            String[] alignmentEntries = StringUtils.split(modelAlignment, ':');
            if (org.apache.commons.lang3.ArrayUtils.getLength(alignmentEntries) == PAIR_SIZE) {
                final String modelName = alignmentEntries[0];
                final List<Fragment> fragmentList = parseFragmentsList(alignmentEntries);
                if (!result.containsKey(modelName)) {
                    result.put(modelName, fragmentList);
                } else {
                    result.get(modelName).addAll(fragmentList);
                }
            }
        }
        return result;
    }

    private static List<Fragment> parseFragmentsList(final String[] alignmentEntries) {
        final List<Fragment> fragmentList = Lists.newArrayList();
        int resNo = 1;
        if (StringUtils.contains(alignmentEntries[1], '|')) {
            final String[] fragments = StringUtils.split(alignmentEntries[1], '|');
            for (String fragment : fragments) {
                final String[] entries = StringUtils.split(fragment, ',');
                fragmentList.add(getFragment(resNo, entries));
                resNo += Integer.parseInt(entries[1]);
            }
        } else {
            final String[] entries = StringUtils.split(alignmentEntries[1], ',');
            fragmentList.add(getFragment(resNo, entries));
            resNo += Integer.parseInt(entries[1]);
        }
        return fragmentList;
    }

    private static Fragment getFragment(final int resNo, final String[] entries) {
        final Fragment frag = new Fragment();
        frag.setResidue(entries[0]);
        frag.setDistance(entries[1]);
        frag.setResNo(String.valueOf(resNo));
        return frag;
    }

}
