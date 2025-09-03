package edu.put.ma.utils;

import java.util.List;

import org.apache.commons.lang3.EnumUtils;

import com.google.common.collect.Lists;

public final class ArrayUtils {

    private ArrayUtils() {
        // hidden constructor
    }

    public static final <E extends Enum<E>> String getEnumNamesString(final Class<E> enumClass) {
        final StringBuilder result = new StringBuilder();
        for (E enumValue : EnumUtils.getEnumList(enumClass)) {
            if (result.length() > 0) {
                result.append(", ");
            }
            result.append(enumValue.toString());
        }
        return result.toString();
    }

    public static final <E extends Enum<E>> List<String> getEnumNames(final Class<E> enumClass) {
        final List<String> result = Lists.newArrayList();
        for (E enumValue : EnumUtils.getEnumList(enumClass)) {
            result.add(enumValue.toString());
        }
        return result;
    }
}
