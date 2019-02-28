Contributing
=======================

We aim to write an efficient and consistent code. So, you should keep in mind
the balance between memory utilization and code efficiency and pay special
attention to cache utilization. This aspect is very important when using
multi-threads applications with shared memory.

Style
-----------------------

Please, follow the general [Julia coding
style](https://docs.julialang.org/en/v1/manual/style-guide). Since it is too
long to describe all details here, study the code already written. However,
in general,

- Name classes, methods, and variables as clear and meaningful as possible;

- Write short commentaries on the code flow to reading more accessible and
  faster;

- Properly document the code, especially the data structures and methods
  definitions. Do not forget to link/refer them;

- 4-space indentation, no trailing spaces, no tabs, Unix/POSIX end of line.
  Try to keep line within 80 columns and do not exceed 90 columns;

- Do not use one-liner branches. Always use `if...end` even it uses
  two more lines. The code must be as clear and easy to read as possible:

```julia
 # Don't do it
a > 1 && b += func(a)

# Ah, way better and clear
if a > 1
    b += func(a)
end
```

- Avoid dense expressions where possible e.g. prefer nested `if`s over complex
  nested `?`s;

- Explicit return should be preferred except in short-form method definitions.
  For functions that do not return values, explicitly use `nothing` at the end;

- Do not use system specific code/headers. Your code must compile in several
  systems with minimum change;

- Make sure that, for each method/function, you write unit tests that cover
  all corner cases, and few regular cases (> 1);

- Do not commit or do a pull request until the code pass in all tests.
