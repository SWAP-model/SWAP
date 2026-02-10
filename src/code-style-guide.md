Here are some guidelines for the coda I researched and I want to stick with.

# Letter Case
It seems like it has been agreed in the community to switch to lowecase everywhere and avoid outdated uppercase. That came from the fact that first computers did not have the lowercase letters, and we are a bit further past that stage.

Add blank lines and whitespace liberally to separate logical sections and improve visual parsing. Always separate terms in arithmetic expressions with spaces (e.g., `z = c*x + d*y`).

# Indentation

Always use spaces, never tabs, for indentation. Tab characters are not part of the Fortran character set and display inconsistently across different editors, terminals, and compilers. Set your editor to emit spaces when you press the tab key.
Indentation Depth

The amount of indentation is a matter of preference, with the most common choices being 2, 3, or 4 spaces per indentation level. Different style guides recommend:

    2 spaces: Recommended by the JULES project and Google's style guides for other languages​

    4 spaces: Used by Fortran stdlib and many contemporary programming languages

we should stick to **4 spaces**.

## What to Indent

Indent all code blocks including:

- Loop bodies (do, do while)
- Conditional blocks (if, else, select case)
- Subroutine and function bodies
- Module and submodule contents

Comments should be indented with the code within their block to maintain visual structure.

There is a discussion on whether to indent the code below `contains`, but the most important is that at least within a module everything remains consistent. For the sake of clarity I would suggest using the indentation to make it easier for folks coming from other languages like Python or JS to read through.

# Naming Conventions
## Choose the right name
Descriptive, meaningful names are essential—avoid cryptic abbreviations. Keep names to one or two syllables when possible. Single-character names should only be used for loop counters or standard mathematical notation (i, j, k).

For mathematical variables, short notation matching literature conventions is acceptable (Ylm, Gamma, Enl). For other entities, prefer clarity over brevity—speed_of_sound(AIR) is better than sos or a.

Use consistent naming patterns throughout your codebase. Common disambiguation strategies include suffixing modules with _m and types with _t, or prefixing library modules with a package identifier.

## Differenciating item groups by case
There are different conventions out there but all agree that most important is to stay consistent. It is great, however, in Python to always be sure that when you see PascalCase it's a class, snake_case it's a function or variable and SCREAMING_SNAKE_CASE is a parameter. Throughout SWAP codebase probably every approach was used.

I will suggest an approach I will try to follow, but it is not really enforcable:

; use snake-case (underscores) for multi-word names (speed_of_sound, temp_wall). This is kind of in line with the fact that Fortran is case insensitive.

# Explicit Declarations

Always use implicit none at the start of every program, module, submodule, and standalone procedure to enforce proper typing (no implicit typing defaults). This prevents implicit typing errors and makes all variable declarations explicit and traceable. Many compilers offer flags like -fimplicit-none to enforce this globally.


# Control Flow and Structure

Avoid goto statements whenever possible. Modern Fortran provides structured alternatives—use named loops with exit and cycle for complex flow control. Use end subroutine foo or end function foo instead of bare end statements to clarify what's being terminated.

Label deeply nested loops so their structure is clear, and prefer one-line if statements for simple conditionals to maintain conciseness.

# Documentation

Include at least one comment line at the beginning of each function or subroutine explaining its purpose. Document procedure arguments either separately or inline with declarations. Capitalize comments as normal prose, not as code.

Ideally though we should use FORD-compliant docstrings. FORD is fortran documentation software (kind of like MKDocs for fortran). This will allow to automatically generate API documentation directly from code.

# Resource Management

Prefer allocatable arrays over pointer arrays when possible, as pointers require explicit deallocation to prevent memory leaks. Use list-directed reads rather than formatted reads for more robust input handling. We will also be using ASSOCIATE to link variable names inside modules with the states that will sit on top of all individual modules.

# functions
I learnt that in a function without specified output variable, the name of the function itself is by default treated as return value. We could follw this for now, because this is how the functions are written for now, but maybe in the future we might reconsider and switch to more explicit style.