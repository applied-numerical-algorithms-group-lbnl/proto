import clang.cindex
import os
import subprocess
from dataclasses import dataclass
from typing import List, Dict, Optional, Set

@dataclass
class FunctionInfo:
    name: str
    signature: str
    return_type: str
    parameters: List[Dict[str, str]]
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
    comment: Optional[str] = None
    is_template: bool = False
    template_parameters: List[str] = None

@dataclass
class NamespaceInfo:
    name: str
    classes: List[ClassInfo]
    functions: List[FunctionInfo]
    nested_namespaces: List['NamespaceInfo']
    source: str

def extract_comment(cursor):
    """Extract documentation comment for a cursor if available."""
    comment = cursor.brief_comment
    if not comment:
        # Try to get raw comment which might include doxygen
        comment = cursor.raw_comment
    return comment

def get_source_text(cursor, tu_source):
    """Extract source code for a given cursor."""
    start = cursor.extent.start
    end = cursor.extent.end
    
    # Get source lines
    start_line = start.line - 1  # 0-based indexing
    end_line = end.line - 1
    
    source_lines = tu_source.splitlines()[start_line:end_line+1]
    
    # Adjust first and last line to respect column offsets
    if source_lines:
        source_lines[0] = source_lines[0][start.column-1:]
        source_lines[-1] = source_lines[-1][:end.column-1]
    
    return '\n'.join(source_lines)

def parse_function(cursor, tu_source):
    """Parse a function declaration."""
    name = cursor.spelling
    signature = cursor.displayname
    return_type = cursor.result_type.spelling
    
    # Extract parameters
    parameters = []
    for param in cursor.get_arguments():
        parameters.append({
            'name': param.spelling,
            'type': param.type.spelling
        })
    
    location = f"{cursor.location.file.name}:{cursor.location.line}"
    source = get_source_text(cursor, tu_source)
    comment = extract_comment(cursor)

    print("Parsed function: " + cursor.spelling)
    return FunctionInfo(
        name=name,
        signature=signature,
        return_type=return_type,
        parameters=parameters,
        location=location,
        source=source,
        comment=comment
    )

def parse_class(cursor, tu_source):
    """Parse a class declaration including its methods and fields."""
    name = cursor.spelling
    location = f"{cursor.location.file.name}:{cursor.location.line}"
    methods = []
    fields = []
    base_classes = []
    is_template = False
    template_parameters = []
    
    # Check if this is a template class
    for child in cursor.get_children():
        if child.kind == clang.cindex.CursorKind.TEMPLATE_TYPE_PARAMETER:
            is_template = True
            template_parameters.append(child.spelling)
    
    # Extract base classes
    for child in cursor.get_children():
        if child.kind == clang.cindex.CursorKind.CXX_BASE_SPECIFIER:
            base_classes.append(child.spelling)
    
    # Extract methods and fields
    for child in cursor.get_children():
        if child.kind in [
            clang.cindex.CursorKind.CXX_METHOD,
            clang.cindex.CursorKind.FUNCTION_TEMPLATE,
            clang.cindex.CursorKind.CONSTRUCTOR,
            clang.cindex.CursorKind.DESTRUCTOR
        ]:
            methods.append(parse_function(child, tu_source))
        elif child.kind == clang.cindex.CursorKind.FIELD_DECL:
            fields.append({
                'name': child.spelling,
                'type': child.type.spelling,
                'comment': extract_comment(child)
            })
    
    source = get_source_text(cursor, tu_source)
    comment = extract_comment(cursor)
    
    print("Parsed class: " + cursor.spelling)
    return ClassInfo(
        name=name,
        location=location,
        source=source,
        methods=methods,
        fields=fields,
        base_classes=base_classes,
        comment=comment,
        is_template=is_template,
        template_parameters=template_parameters
    )


def parse_namespace(cursor, tu_source):
    """Parse a namespace and its contents recursively."""
    name = cursor.spelling or "<anonymous>"
    classes = []
    functions = []
    nested_namespaces = []
    
    for child in cursor.get_children():
        if child.kind == clang.cindex.CursorKind.NAMESPACE:
            nested_namespaces.append(parse_namespace(child, tu_source))
        elif child.kind in [clang.cindex.CursorKind.CLASS_DECL, clang.cindex.CursorKind.STRUCT_DECL, clang.cindex.CursorKind.CLASS_TEMPLATE]:
            classes.append(parse_class(child, tu_source))
        elif child.kind in [clang.cindex.CursorKind.FUNCTION_DECL, clang.cindex.CursorKind.FUNCTION_TEMPLATE]:
            functions.append(parse_function(child, tu_source))
    
    source = get_source_text(cursor, tu_source)
    print("Parsed namespace: " + cursor.spelling)
    return NamespaceInfo(
        name=name,
        classes=classes,
        functions=functions,
        nested_namespaces=nested_namespaces,
        source=source
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

def parse_header_file(file_path, include_paths=None, compiler_name='clang++'):
    """Parse a C++ header file and extract its structure."""
    index = clang.cindex.Index.create()
    args = setup_compiler_args(include_paths, compiler_name)

    # Parse the file
    try:
        translation_unit = index.parse(file_path, args=args)
        
        # Check for parsing errors
        for diag in translation_unit.diagnostics:
            if diag.severity >= clang.cindex.Diagnostic.Error:
                print(f"Error parsing {file_path}: {diag.spelling}")
        
        # Get file content as string for source extraction
        with open(file_path, 'r') as f:
            file_content = f.read()
        
        # Root cursor represents the translation unit
        root = translation_unit.cursor
        
        # Extract top-level declarations
        global_namespace = NamespaceInfo(
            name="",
            classes=[],
            functions=[],
            nested_namespaces=[],
            source=""
        )
        
        for child in root.get_children():
            # Filter out declarations from other files
            if child.location.file and child.location.file.name == file_path:
                if child.kind == clang.cindex.CursorKind.NAMESPACE:
                    global_namespace.nested_namespaces.append(parse_namespace(child, file_content))
                elif child.kind in [clang.cindex.CursorKind.CLASS_DECL, clang.cindex.CursorKind.STRUCT_DECL, clang.cindex.CursorKind.CLASS_TEMPLATE]:
                    global_namespace.classes.append(parse_class(child, file_content))
                elif child.kind in [clang.cindex.CursorKind.FUNCTION_DECL, clang.cindex.CursorKind.FUNCTION_TEMPLATE]:
                    global_namespace.functions.append(parse_function(child, file_content))
        
        return global_namespace
        
    except Exception as e:
        print(f"Error parsing {file_path}: {str(e)}")
        return None



def parse_codebase(root_dir, include_paths=None):
    """Parse all header files in a directory recursively."""
    all_parsed_files = {}
    main_filepath = os.path.join(root_dir, 'Main.H')
    all_parsed_files[main_filepath] = parse_header_file(main_filepath, include_paths)
    # for dirpath, _, filenames in os.walk(root_dir):
    #     for filename in filenames:
    #         if filename.endswith(('.h', '.hpp', '.hxx', '.H')):
    #             file_path = os.path.join(dirpath, filename)
    #             parsed_file = parse_header_file(file_path, include_paths)
    #             if parsed_file:
    #                 all_parsed_files[file_path] = parsed_file
    
    return all_parsed_files

if __name__ == "__main__":
    include_paths = [
        "./include/base",
        "./include/base/implem"
    ]
    codebase_structure = parse_codebase("./include")

    #Print summary of parsed structures
    for file_path, structure in codebase_structure.items():
        print(f"File: {file_path}")
        
        # Count classes, functions, etc.
        classes_count = len(structure.classes)
        for ns in structure.nested_namespaces:
            classes_count += len(ns.classes)
            
        functions_count = len(structure.functions)
        for ns in structure.nested_namespaces:
            functions_count += len(ns.functions)
            
        print(f"  Classes: {classes_count}")
        print(f"  Functions: {functions_count}")
        print()