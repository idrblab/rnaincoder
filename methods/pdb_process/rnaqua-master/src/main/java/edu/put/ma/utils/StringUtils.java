package edu.put.ma.utils;

public final class StringUtils {

    public static final char NEW_LINE = '\n';

    public static final boolean ALIGNED_LEFT_MODE = true;

    private StringUtils() {
        // hidden constructor
    }

    public static final String deleteLastCharacter(final String input) {
        final int inputLength = org.apache.commons.lang3.StringUtils.length(input);
        if (inputLength > 1) {
            return org.apache.commons.lang3.StringUtils.substring(input, 0, inputLength - 1);
        }
        return input;
    }

    public static final String extendString(final String input, final int expectedSize, final char character,
            final boolean alignedLeft) {
        final int inputLength = org.apache.commons.lang3.StringUtils.length(input);
        final StringBuilder result = new StringBuilder();
        if (alignedLeft) {
            result.append(input).append(
                    org.apache.commons.lang3.StringUtils.repeat(character, expectedSize - inputLength));
        } else {
            result.append(org.apache.commons.lang3.StringUtils.repeat(character, expectedSize - inputLength))
                    .append(input);
        }
        return result.toString();
    }

    public static final String deleteDanglingChars(final String input, final String character) {
        String output = org.apache.commons.lang3.StringUtils.trim(input);
        if (org.apache.commons.lang3.StringUtils.isNotBlank(output)) {
            if (output.endsWith(character)) {
                output = output.substring(0, output.length() - 1);
            }
            if (output.startsWith(character)) {
                output = output.substring(1);
            }
        }
        return output;
    }

    public static final void extendStringBuilder(final StringBuilder sb, final String value,
            final String label) {
        if ((sb != null) && (org.apache.commons.lang3.StringUtils.isNotBlank(value))) {
            sb.append(String.format("\n%s: ", label)).append(value);
        }
    }

    public static final String deleteAllSpaces(final String inputString) {
        if (org.apache.commons.lang3.StringUtils.isNotBlank(inputString)) {
            return inputString.replaceAll(" +", "");
        }
        return inputString;
    }
}
