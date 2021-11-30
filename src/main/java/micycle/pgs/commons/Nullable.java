package micycle.pgs.commons;

import java.lang.annotation.Documented;
import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * An element annotated with Nullable claims <code>null</code> value is
 * perfectly <em>valid</em> to return (for methods), pass to (parameters) and
 * hold (local variables and fields).
 */
@Documented
@Retention(RetentionPolicy.CLASS)
@Target({ ElementType.METHOD, ElementType.FIELD, ElementType.PARAMETER, ElementType.LOCAL_VARIABLE })
public @interface Nullable {
	String value() default "";
}