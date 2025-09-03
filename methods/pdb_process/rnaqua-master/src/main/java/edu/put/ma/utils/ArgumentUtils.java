package edu.put.ma.utils;

import java.util.Arrays;

import org.apache.commons.lang3.ArrayUtils;

import com.google.common.collect.ImmutableSet;

public final class ArgumentUtils {

    private static final int SINGLE_ARGUMENT_SIZE = 2;

    private ArgumentUtils() {
        // hidden constructor
    }

    public static final String[] retrieveArgByOpt(final String[] args, final String opt) {
        final int argIndex = getIndexOfArgByOpt(args, opt);
        if (argIndex >= 0) {
            return Arrays.copyOfRange(args, argIndex, argIndex + SINGLE_ARGUMENT_SIZE);
        }
        return ArrayUtils.EMPTY_STRING_ARRAY;
    }

    public static final String[] retrieveArgByOpt(final String[] args, final ImmutableSet<String> opt) {
        String[] result = ArrayUtils.EMPTY_STRING_ARRAY;
        for (String optCode : opt) {
            final String[] array = retrieveArgByOpt(args, optCode);
            if (array != result) {
                result = array;
                break;
            }
        }
        return result;
    }

    public static final String[] removeArgByOpt(final String[] args, final String opt) {
        final int argsCount = ArrayUtils.getLength(args);
        final int argIndex = getIndexOfArgByOpt(args, opt);
        if ((argIndex >= 0) && (argIndex < argsCount)) {
            String[] prefix = null;
            if (argIndex > 0) {
                prefix = Arrays.copyOfRange(args, 0, argIndex);
            }
            String[] postfix = null;
            if (argIndex + SINGLE_ARGUMENT_SIZE < argsCount) {
                postfix = Arrays.copyOfRange(args, argIndex + SINGLE_ARGUMENT_SIZE, argsCount);
            }
            return ArrayUtils.addAll(prefix, postfix);
        }
        return ArrayUtils.EMPTY_STRING_ARRAY;
    }

    public static final String[] removeArgByOpt(final String[] args, final ImmutableSet<String> opt) {
        String[] result = ArrayUtils.EMPTY_STRING_ARRAY;
        for (String optCode : opt) {
            final String[] array = removeArgByOpt(args, optCode);
            if (array != result) {
                result = array;
                break;
            }
        }
        return result;
    }

    private static final int getIndexOfArgByOpt(final String[] args, final String opt) {
        int result = -1;
        final String shortOptionCode = "-" + opt;
        final String longOptionCode = "--" + opt;
        PreconditionUtils.checkIfStringIsBlank(opt, "Option code");
        final int argsCount = ArrayUtils.getLength(args);
        for (int argIndex = 0; argIndex < argsCount; argIndex++) {
            if ((shortOptionCode.equalsIgnoreCase(args[argIndex]))
                    || (longOptionCode.equalsIgnoreCase(args[argIndex]))) {
                result = argIndex;
                break;
            }
        }
        return result;
    }
}
