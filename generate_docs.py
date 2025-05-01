import sys
import traceback
import clang.cindex
import os
import subprocess
from dataclasses import dataclass
from typing import Any, List, Dict, Optional

# =================================================================================================
# CODE PARSING 
@dataclass
class FunctionInfo:
    name: str
    signature: str
    return_type: str
    parameters: List[Dict[str, Any]]
    location: str
    attributes: Dict[str,str]
    source: str
    comment: Optional[str] = None
@dataclass
class EnumInfo:
    name: str
    values: List[Dict[str, str]]
    location: str
    source: str
    comment: Optional[str] = None
@dataclass
class TypedefInfo:
    name: str
    base_type: str
    location: str
    source: str
    comment: Optional[str] = None
@dataclass
class ClassInfo:
    name: str
    location: str
    source: str
    methods: List[FunctionInfo]
    fields: List[Dict[str, str]]
    base_classes: List[str]
    friend_functions: List[FunctionInfo] = None
    friend_classes: List[str] = None
    comment: Optional[str] = None
    is_template: bool = False
    template_parameters: List[str] = None

@dataclass
class NamespaceInfo:
    name: str
    classes: List[ClassInfo]
    functions: List[FunctionInfo]
    enums: List[EnumInfo]
    typedefs: List[TypedefInfo]
    nested_namespaces: List['NamespaceInfo']
    source: str
    using_directives: List[str] = None

def extract_comment(cursor):
    comment = cursor.brief_comment
    if not comment:
        # Try to get raw comment which might include doxygen
        comment = cursor.raw_comment
    return comment

def get_source_text(cursor, tu_source):
    start = cursor.extent.start
    end = cursor.extent.end
    
    start_line = start.line - 1  # 0-based indexing
    end_line = end.line - 1
    
    source_lines = tu_source.splitlines()[start_line:end_line+1]
    
    if source_lines:
        source_lines[0] = source_lines[0][start.column-1:]
        source_lines[-1] = source_lines[-1][:end.column-1]
    
    return '\n'.join(source_lines)

def is_method(cursor):
    if cursor.kind in [
        clang.cindex.CursorKind.CXX_METHOD,
        clang.cindex.CursorKind.FUNCTION_TEMPLATE,
        clang.cindex.CursorKind.CONSTRUCTOR,
        clang.cindex.CursorKind.DESTRUCTOR
    ]:
        return True
    return False

def parse_function(cursor, tu_source):
    name = cursor.spelling
    signature = cursor.displayname
    return_type = cursor.result_type.spelling

# Safe attribute extraction with fallbacks
    attributes = {}
    
    # Storage class (static, extern, etc.)
    try:
        attributes['storage_class'] = cursor.storage_class.name
    except (AttributeError, ValueError):
        attributes['storage_class'] = None
    
    # Check if method is inline
    try:
        attributes['is_inline'] = cursor.is_inline()
    except (AttributeError, ValueError):
        # Fallback: check tokens for 'inline' keyword
        attributes['is_inline'] = 'inline' in get_source_text(cursor, tu_source).split()
    
    # For methods only (not free functions)
    if cursor.kind in [
        clang.cindex.CursorKind.CXX_METHOD,
        clang.cindex.CursorKind.CONSTRUCTOR,
        clang.cindex.CursorKind.DESTRUCTOR
    ]:
        # Method-specific attributes
        try:
            attributes['is_const'] = cursor.is_const_method()
        except (AttributeError, ValueError):
            attributes['is_const'] = False
            
        try:
            attributes['is_virtual'] = cursor.is_virtual_method()
        except (AttributeError, ValueError):
            attributes['is_virtual'] = False
            
        try:
            attributes['is_pure_virtual'] = cursor.is_pure_virtual_method()
        except (AttributeError, ValueError):
            attributes['is_pure_virtual'] = False
            
        try:
            attributes['is_static'] = cursor.is_static_method()
        except (AttributeError, ValueError):
            attributes['is_static'] = False
            
        # Access specifier
        try:
            attributes['access_specifier'] = cursor.access_specifier.name.lower()
        except (AttributeError, ValueError):
            attributes['access_specifier'] = None
    
    # For all function types
    # Check for noexcept (libclang doesn't have direct API)
    source_text = get_source_text(cursor, tu_source)
    attributes['is_noexcept'] = 'noexcept' in source_text
    attributes['is_constexpr'] = 'constexpr' in source_text.split()
    
    # Look for override, final keywords
    attributes['is_override'] = 'override' in source_text.split()
    attributes['is_final'] = 'final' in source_text.split()
    # Extract parameters
    parameters = []
    for param in cursor.get_arguments():
        default_value = None
        for child in param.get_children():
            # Look for default value expression
            if child.kind in [clang.cindex.CursorKind.INTEGER_LITERAL, 
                             clang.cindex.CursorKind.FLOATING_LITERAL,
                             clang.cindex.CursorKind.STRING_LITERAL,
                             clang.cindex.CursorKind.CHARACTER_LITERAL,
                             clang.cindex.CursorKind.CXX_BOOL_LITERAL_EXPR,
                             clang.cindex.CursorKind.CXX_NULL_PTR_LITERAL_EXPR,
                             clang.cindex.CursorKind.UNEXPOSED_EXPR]:
                default_value = get_source_text(child, tu_source)
                break
        
        parameters.append({
            'name': param.spelling,
            'type': param.type.spelling,
            'default_value': default_value
        })
    
    location = f"{cursor.location.file.name}:{cursor.location.line}"
    source = get_source_text(cursor, tu_source)
    comment = extract_comment(cursor)

    return FunctionInfo(
        name=name,
        signature=signature,
        return_type=return_type,
        attributes=attributes,
        parameters=parameters,
        location=location,
        source=source,
        comment=comment
    )

def parse_class(cursor, tu_source):
    name = cursor.spelling
    location = f"{cursor.location.file.name}:{cursor.location.line}"
    methods = []
    fields = []
    base_classes = []
    is_template = False
    template_parameters = []
    friend_functions = []
    friend_classes=[]
    for child in cursor.get_children():
        if child.kind == clang.cindex.CursorKind.TEMPLATE_TYPE_PARAMETER:
            is_template = True
            template_parameters.append(child.spelling)
        elif child.kind == clang.cindex.CursorKind.CXX_BASE_SPECIFIER:
            base_classes.append(child.spelling)
        elif is_method(child):
            methods.append(parse_function(child, tu_source))
        elif child.kind == clang.cindex.CursorKind.FIELD_DECL:
            fields.append({
                'name': child.spelling,
                'type': child.type.spelling,
                'comment': extract_comment(child)
            })
        elif child.kind == clang.cindex.CursorKind.FRIEND_DECL:
            for grandchild in child.get_children():
                if grandchild.kind in [clang.cindex.CursorKind.FUNCTION_DECL, clang.cindex.CursorKind.FUNCTION_TEMPLATE]:
                    friend_functions.append(parse_function(grandchild, tu_source))
                elif grandchild.kind in [clang.cindex.CursorKind.CLASS_DECL, clang.cindex.CursorKind.STRUCT_DECL]:
                    friend_classes.append(grandchild.spelling)
    
    source = get_source_text(cursor, tu_source)
    comment = extract_comment(cursor)
    
    return ClassInfo(
        name=name,
        location=location,
        source=source,
        methods=methods,
        fields=fields,
        base_classes=base_classes,
        friend_functions=friend_functions,
        friend_classes=friend_classes,
        comment=comment,
        is_template=is_template,
        template_parameters=template_parameters
    )

def parse_enum(cursor, tu_source):
    name = cursor.spelling
    location = f"{cursor.location.file.name}:{cursor.location.line}"
    values = []
    
    for child in cursor.get_children():
        if child.kind == clang.cindex.CursorKind.ENUM_CONSTANT_DECL:
            # Try to get the enum value
            value = None
            for grandchild in child.get_children():
                if grandchild.kind == clang.cindex.CursorKind.INTEGER_LITERAL:
                    value = grandchild.spelling
                    break
            
            values.append({
                'name': child.spelling,
                'value': value,
                'comment': extract_comment(child)
            })
    
    source = get_source_text(cursor, tu_source)
    comment = extract_comment(cursor)
    
    return EnumInfo(
        name=name,
        location=location,
        source=source,
        values=values,
        comment=comment
    )

def parse_typedef(cursor, tu_source):
    name = cursor.spelling
    base_type = cursor.underlying_typedef_type.spelling
    location = f"{cursor.location.file.name}:{cursor.location.line}"
    source = get_source_text(cursor, tu_source)
    comment = extract_comment(cursor)
    
    return TypedefInfo(
        name=name,
        base_type=base_type,
        source=source,
        location=location,
        comment=comment
    )

def parse_namespace(cursor, tu_source):
    name = cursor.spelling or "<anonymous>"
    classes = []
    functions = []
    nested_namespaces = []
    enums = []
    typedefs = []
    using_directives = []
    for child in cursor.get_children():
        if child.kind == clang.cindex.CursorKind.NAMESPACE:
            nested_namespaces.append(parse_namespace(child, tu_source))
        elif child.kind in [clang.cindex.CursorKind.CLASS_DECL, clang.cindex.CursorKind.STRUCT_DECL, clang.cindex.CursorKind.CLASS_TEMPLATE]:
            classes.append(parse_class(child, tu_source))
        elif child.kind in [clang.cindex.CursorKind.FUNCTION_DECL, clang.cindex.CursorKind.FUNCTION_TEMPLATE]:
            functions.append(parse_function(child, tu_source))
        elif child.kind == clang.cindex.CursorKind.ENUM_DECL:
            enums.append(parse_enum(child, tu_source))
        elif child.kind == clang.cindex.CursorKind.TYPEDEF_DECL:
            typedefs.append(parse_typedef(child, tu_source))
        elif child.kind == clang.cindex.CursorKind.USING_DIRECTIVE:
            using_directives.append(child.spelling)
    
    source = get_source_text(cursor, tu_source)
    return NamespaceInfo(
        name=name,
        classes=classes,
        functions=functions,
        enums=enums,
        typedefs=typedefs,
        nested_namespaces=nested_namespaces,
        source=source,
        using_directives=using_directives
    )

def get_system_includes(compiler_name='clang++'):
    
    try:
        # Run preprocesser to extract stderr
        cmd = [compiler_name, '-E', '-x', 'c++', '-v', '/dev/null']
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        begin_includes_marker = "#include <...> search starts here:"
        end_includes_marker = "End of search list."

        # Extract paths between specific markers in the output
        includes = []
        in_includes = False
        for line in result.stderr.splitlines():
            if line.strip() == begin_includes_marker:
                in_includes = True
                continue
            elif line.strip() == end_includes_marker:
                in_includes = False
                continue
            
            if in_includes:
                includes.append(line.strip())
                
        return includes
    except Exception as e:
        print(f"Warning: Could not detect system includes: {e}")
        return []

def setup_compiler_args(include_paths=None, compiler_name='clang++'):
    args = [
        '-x', 'c++',
        '-std=c++17', 
        '-Wno-deprecated-declarations',
        '-ferror-limit=0',
    ]
    
    # Add system includes
    system_includes = get_system_includes(compiler_name)
    for path in system_includes:
        args.extend(['-isystem', path])
    
    # Add user include paths
    if include_paths:
        for path in include_paths:
            args.append(f'-I{path}')

    return args

def get_method_source(cursor, file_content):
    """Extract the source code for a method from file content."""
    start = cursor.extent.start.offset
    end = cursor.extent.end.offset
    
    # Make sure offsets are valid
    if start >= 0 and end >= 0 and start < len(file_content) and end <= len(file_content):
        return file_content[start:end]
    return ""

def find_implementations(cursor, implem_file_path, class_data=None, func_data=None):
    implems = {}
    with open(implem_file_path, 'r') as f:
        implem_file_content = f.read()

    for child in cursor.get_children():
        if not child.location.file:
            continue
        
        if child.location.file.name == implem_file_path:
            if class_data and child.kind in [
                clang.cindex.CursorKind.CXX_METHOD,
                clang.cindex.CursorKind.FUNCTION_TEMPLATE,
                clang.cindex.CursorKind.CONSTRUCTOR,
                clang.cindex.CursorKind.DESTRUCTOR
            ]:
                if child.semantic_parent and child.semantic_parent.spelling == class_data.name:
                    impl = get_method_source(child, implem_file_content)
                    implems[child.spelling] = impl
            
            elif func_data and child.spelling == func_data.name:
                impl = get_method_source(child, implem_file_content)
                implems[child.spelling] = impl

            elif child.kind == clang.cindex.CursorKind.NAMESPACE:
                nested_implems = find_implementations(child, implem_file_path, class_data, func_data)
                implems.update(nested_implems)

    if class_data:
        for method_data in class_data.methods:
            if method_data.name in implems:
                method_data.source = implems[method_data.name]
    if func_data:
        func_data.source = implems[func_data.name]
    return implems

def parse_header_file(file_path, implem_file_path, include_paths=None, compiler_name='clang++'):
    index = clang.cindex.Index.create()
    args = setup_compiler_args(include_paths, compiler_name)

    # Parse the file
    try:
        header_tu = index.parse(file_path, args=args)

        # Check for parsing errors
        for diag in header_tu.diagnostics:
            if diag.severity >= clang.cindex.Diagnostic.Error:
                message = f"Error: {diag.spelling}"
            
                # Add location information if available
                if diag.location.file:
                    message += f" at {diag.location.file.name}:{diag.location.line}:{diag.location.column}"
                
                # Add any additional notes/fixits
                if diag.fixits:
                    for fixit in diag.fixits:
                        message += f"\n  Suggested fix: {fixit.value}"
                
                # Print all diagnostics, but flag errors
                print(message)
        
        # Get file content as string for source extraction
        with open(file_path, 'r') as f:
            header_file_content = f.read()
        
        # Root cursor represents the translation unit
        root = header_tu.cursor
        
        # Extract top-level declarations
        global_namespace = NamespaceInfo(
            name="",
            classes=[],
            functions=[],
            nested_namespaces=[],
            enums=[],
            typedefs=[],
            source=""
        )
        
        for child in root.get_children():
            # Filter out declarations from other files
            if child.location.file and child.location.file.name == file_path:
                if child.kind == clang.cindex.CursorKind.NAMESPACE:
                    global_namespace.nested_namespaces.append(parse_namespace(child, header_file_content))
                elif child.kind in [clang.cindex.CursorKind.CLASS_DECL, clang.cindex.CursorKind.STRUCT_DECL, clang.cindex.CursorKind.CLASS_TEMPLATE]:
                    global_namespace.classes.append(parse_class(child, header_file_content))
                elif child.kind in [clang.cindex.CursorKind.FUNCTION_DECL, clang.cindex.CursorKind.FUNCTION_TEMPLATE]:
                    global_namespace.functions.append(parse_function(child, header_file_content))
                elif child.kind == clang.cindex.CursorKind.ENUM_DECL:
                    global_namespace.enums.append(parse_enum(child, header_file_content))
                elif child.kind == clang.cindex.CursorKind.TYPEDEF_DECL:
                    global_namespace.typedefs.append(parse_typedef(child, header_file_content))
        
        #get implementations
        if implem_file_path and os.path.exists(implem_file_path):
            implem_tu = index.parse(implem_file_path, args=args)
            
            for diag in implem_tu.diagnostics:
                if diag.severity >= clang.cindex.Diagnostic.Error:
                    print(f"Error parsing {implem_file_path}: {diag.spelling}")

            for ci in global_namespace.classes:
                implems = find_implementations(implem_tu.cursor, implem_file_path, ci, None)
            for fi in global_namespace.functions:
                implems = find_implementations(implem_tu.cursor, implem_file_path, None, fi)
            for ns in global_namespace.nested_namespaces:
                for ci in ns.classes:
                    implems = find_implementations(implem_tu.cursor, implem_file_path, ci, None)
                for fi in ns.functions:
                    implems = find_implementations(implem_tu.cursor, implem_file_path, None, fi)

        return global_namespace
        
    except Exception as e:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        print(f"\nError parsing {file_path}:")
        print(f"  Exception type: {exc_type.__name__}")
        print(f"  Exception message: {str(e)}")
        return None

def parse_codebase(root_dir, include_paths=None):
    all_parsed_files = {}
    main_filepath = os.path.join(root_dir, 'Main.H')
    all_parsed_files[main_filepath] = parse_header_file(main_filepath, None, include_paths)
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.endswith(('.h', '.hpp', '.hxx', '.H')):
                file_path = os.path.join(dirpath, filename)
                name, ext = os.path.splitext(filename)
                implem_filename = f"{name}Implem{ext}"
                implem_file_path = f"{dirpath}/implem/{implem_filename}"
                parsed_file = parse_header_file(file_path, implem_file_path, include_paths)
                if parsed_file:
                    all_parsed_files[file_path] = parsed_file
                else:
                    print(f"Warning: Skipping {file_path} due to parsing errors")
    
    return all_parsed_files

def print_structure(codebase_structure):
    for file_path, structure in codebase_structure.items():
        print(f"source file: {file_path}")
        for ns in structure.nested_namespaces:
            for td in ns.typedefs:
                print(f"\tTypedef {td.source}")
            for ei in ns.enums:
                print(f"\tEnum {ei.name}")
                for ev in ei.values:
                    print(f"\t\t{ev['name']}")
            for ci in ns.classes:
                print(f"\tClass {ci.name}")
                print(ci.source)
                for vi in ci.fields:
                    print(f"\t\tField {vi['name']}")
                for fi in ci.methods:
                    print(f"\t\tMethod {fi.name}")
                    # print(fi.source)

# =================================================================================================
# DOCUMENT GENERATION

def create_class_documentation_prompt(class_info, namespace_context=""):

    # Create fully qualified name
    qualified_name = f"{namespace_context}::{class_info.name}" if namespace_context else class_info.name
    
    prompt = f"""
# C++ Class Documentation Request

Generate comprehensive markdown documentation for the following C++ class:

## Basic Information
- **Class Name**: {class_info.name}
- **Fully Qualified Name**: {qualified_name}
- **File Location**: {class_info.location}

## Class Definition
```cpp
{class_info.source}
```
"""
    # Add template information if present
    if class_info.is_template and class_info.template_parameters:
        prompt += "\n## Template Parameters\n"
        for param in class_info.template_parameters:
            prompt += f"- {param}\n"

    # Add inheritance information
    if class_info.base_classes:
        prompt += "\n## Inheritance\n"
        for base in class_info.base_classes:
            prompt += f"- Inherits from: {base}\n"

    # Include existing documentation comment if available
    if class_info.comment:
        prompt += f"\n## Existing Documentation\n{class_info.comment}\n"

    # Add field information
    if class_info.fields:
        prompt += "\n## Member Variables\n"
        for field in class_info.fields:
            field_info = f"- **{field['name']}**: {field['type']}"
            if field.get('comment'):
                field_info += f"\n  Comment: {field['comment']}"
            prompt += field_info + "\n"

    # Add method signatures (detailed method docs will be separate)
    if class_info.methods:
        prompt += "\n## Methods Overview\n"
        for method in class_info.methods:
            prompt += f"- {method.signature}\n"
    
    # Instructions for the LLM
# Instructions for the LLM
    prompt += """
## Documentation Requirements

Please provide clear, comprehensive documentation for this class that includes:

1. **Class Purpose and Overview**: A concise description of what this class represents and its role.

2. **Usage Examples**: At least one example showing typical usage of this class.

3. **Design Patterns**: Identify any design patterns or principles employed.

4. **Thread Safety**: Mention any thread safety considerations if applicable.

5. **Member Variables Documentation**: For each member variable, explain:
   - Its purpose and role in the class
   - Any constraints or valid values
   - Its relationship to class functionality

6. **Methods Organization**: Group related methods together under logical sections.

7. **Key Method Explanations**: For important methods, provide:
   - Purpose and functionality
   - How they should be used
   - Any preconditions and postconditions
   - Error handling behavior
   
8. **Best Practices**: Any recommendations for effectively using this class.

Format the documentation using Markdown with appropriate headers, code blocks, and emphasis.
"""
    return prompt

if __name__ == "__main__":
    include_paths = [
        "./include/base",
        "./include/base/implem"
    ]
    codebase_structure = parse_codebase("./include", include_paths)
    print_structure(codebase_structure)
    for file_path, structure in codebase_structure.items():
        context = ""
        for ns in structure.nested_namespaces:
            context += f"::{ns.name}"
            for ci in ns.classes:
                prompt = create_class_documentation_prompt(ci, context)
                print(prompt)
    