using Printf

function parse_namelist(filename::AbstractString)
    # read file, strip comments and empty lines
    lines = String[]
    open(filename) do f
        for line in eachline(f)
            # remove anything after a comment character
            if occursin("!", line)
                line = split(line, "!")[1]
            end
            line = strip(line)
            if !isempty(line)
                push!(lines, line)
            end
        end
    end

    groups = Dict{String, Dict{String, Any}}()
    i = 1
    n = length(lines)

    while i <= n
        line = lines[i]
        if startswith(line, "&")
            # start of a namelist group
            group_name = strip(line[2:end])
            group_lines = String[]
            i += 1
            while i <= n
                l = lines[i]
                if occursin("/", l)
                    # split at '/', keep the part before it
                    parts = split(l, "/")
                    if !isempty(parts[1])
                        push!(group_lines, parts[1])
                    end
                    break
                else
                    push!(group_lines, l)
                    i += 1
                end
            end
            # parse the content of this group
            groups[group_name] = parse_group(group_lines)
        end
        i += 1
    end
    return groups
end

function parse_group(lines::Vector{String})
    # join all lines with spaces to handle multi‑line values
    text = join(lines, " ")
    tokens = tokenize(text)
    dict = Dict{String, Any}()

    i = 1
    while i <= length(tokens)
        tok = tokens[i]
        if tok isa String && tok != "="
            # variable name
            name = tok
            i += 1
            if i > length(tokens) || tokens[i] != "="
                error("Expected '=' after variable name '$name'")
            end
            i += 1
            # collect value tokens until the next variable name
            val_tokens = String[]
            while i <= length(tokens)
                t = tokens[i]
                # heuristic: a token that looks like a variable name and
                # is not a known boolean literal -> start of next assignment
                if is_variable_name(t) && !is_boolean_literal(t)
                    break
                else
                    push!(val_tokens, t)
                    i += 1
                end
            end
            dict[name] = parse_value(val_tokens)
        else
            i += 1
        end
    end
    return dict
end

function tokenize(text::String)
    tokens = Vector{String}()
    i = 1
    n = length(text)
    while i <= n
        c = text[i]
        if isspace(c)
            i += 1
            continue
        elseif c == '=' || c == ',' || c == '/' || c == '&'
            push!(tokens, string(c))
            i += 1
        elseif c == '\'' || c == '"'
            # quoted string
            quote_char = c
            i += 1
            start = i
            while i <= n && text[i] != quote_char
                # skip escaped quotes (optional: handle \' etc.)
                i += 1
            end
            if i > n
                error("Unclosed string literal")
            end
            str = text[start:i-1]
            push!(tokens, str)      # store without quotes
            i += 1                  # skip closing quote
        elseif c == '.' || c == '-' || c == '+' || isdigit(c)
            # number (integer/float) or logical (e.g. .true.)
            start = i
            while i <= n && (isdigit(text[i]) || text[i] == '.' ||
                             text[i] == 'e' || text[i] == 'E' ||
                             text[i] == '-' || text[i] == '+')
                i += 1
            end
            token = text[start:i-1]
            push!(tokens, token)
        elseif isletter(c) || c == '_'
            # identifier (variable name)
            start = i
            while i <= n && (isletter(text[i]) || isdigit(text[i]) || text[i] == '_')
                i += 1
            end
            token = text[start:i-1]
            push!(tokens, token)
        else
            # skip any unknown characters
            i += 1
        end
    end
    return tokens
end

function parse_value(tokens::Vector{String})
    # if any comma is present, it is an array
    if any(t -> t == ",", tokens)
        values = []
        current = String[]
        for tok in tokens
            if tok == ","
                if !isempty(current)
                    push!(values, parse_scalar(join(current, "")))
                    current = []
                end
            else
                push!(current, tok)
            end
        end
        if !isempty(current)
            push!(values, parse_scalar(join(current, "")))
        end
        return values
    else
        return parse_scalar(join(tokens, ""))
    end
end

function parse_scalar(str::String)
    s = strip(str)
    if isempty(s)
        return nothing
    end
    # boolean
    if s in (".true.", ".TRUE.", "true", "TRUE", "T", "t")
        return true
    elseif s in (".false.", ".FALSE.", "false", "FALSE", "F", "f")
        return false
    end
    # integer
    if tryparse(Int, s) !== nothing
        return parse(Int, s)
    end
    # float (also catches numbers like 5.62000)
    if tryparse(Float64, s) !== nothing
        return parse(Float64, s)
    end
    # otherwise it is a plain string (already unquoted)
    return s
end

# helper functions
is_variable_name(s::String) = occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", s)
is_boolean_literal(s::String) = s in (".true.", ".TRUE.", "true", "TRUE", "T", "t",
                                      ".false.", ".FALSE.", "false", "FALSE", "F", "f")

#groups = parse_namelist("inp01")
