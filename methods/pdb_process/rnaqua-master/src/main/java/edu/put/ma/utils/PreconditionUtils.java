package edu.put.ma.utils;

import java.io.File;
import java.util.List;

import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.lang3.StringUtils;

public final class PreconditionUtils {

    private PreconditionUtils() {
        // hidden constructor
    }

    public static final void checkIfIndexInRange(final int index, final int min, final int max,
            final String prefix) {
        if (!((index >= min) && (index < max))) {
            throw new IndexOutOfBoundsException(String.format("%s index '%s' out of range [%s:%s]", prefix,
                    index, min, max));
        }
    }

    public static void checkIfFileExistsAndIsNotADirectory(final File file, final String prefix) {
        if (!((file.exists()) && (!file.isDirectory()))) {
            throw new IllegalArgumentException(String.format("%s ('%s') doesn't exist or is a directory",
                    prefix, file.getAbsolutePath()));
        }
    }

    public static void checkIfDirectoryExists(final File directory, final String prefix) {
        if (!((directory.exists()) && (directory.isDirectory()))) {
            throw new IllegalArgumentException(String.format("%s ('%s') doesn't exist or isn't a directory",
                    prefix, directory.getAbsolutePath()));
        }
    }

    public static void checkIfFileDoesNotExistOrIsNotADirectory(final File file, final String prefix) {
        if ((file.exists()) && (file.isDirectory())) {
            throw new IllegalArgumentException(String.format("%s ('%s') cannot indicate on a directory", prefix,
                    file.getAbsolutePath()));
        }
    }

    public static final void checkIfStringIsBlank(final String value, final String prefix) {
        if (StringUtils.isBlank(value)) {
            throw new IllegalArgumentException(String.format("%s cannot be blank and should be defined",
                    prefix));
        }
    }

    public static <T> void checkIfListIsEmpty(final List<T> list, final String prefix) {
        if (CollectionUtils.sizeIsEmpty(list)) {
            throw new IllegalArgumentException(String.format("%s cannot be empty", prefix));
        }
    }
}
