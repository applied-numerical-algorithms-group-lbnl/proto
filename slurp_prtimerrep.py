import re
from typing import Dict, List, Any, Optional

def parse_timing_report(content: str) -> Dict[str, Any]:
    """
    Parse the hierarchical timing report into a structured dictionary.
    Stops at detailed breakdown section.
    """
    lines = content.strip().split('\n')

    # Find the start of the actual timing data
    start_idx = 0
    for i, line in enumerate(lines):
        if line.strip().startswith('[0]'):
            start_idx = i
            break
    
    # Find the end of the hierarchical section
    end_idx = len(lines)
    for i in range(start_idx, len(lines)):
        line = lines[i].strip()
        
        # Stop at separator lines or detailed breakdown sections
        if (line.startswith('=======') or 
            line.startswith('-------') or 
            line.startswith('[0]root') or  # Start of detailed breakdown
            (line and not line.startswith('[') and not line.startswith(' '))):
            
            # Check if this is actually a separator (not just a line without brackets)
            if ('=' in line and len(line) > 10) or ('-' in line and len(line) > 10):
                end_idx = i
                break
            # Also stop if we see the detailed breakdown format
            elif line.startswith('[0]root') and i > start_idx + 5:  # Make sure it's not the first occurrence
                end_idx = i
                break
    
    # Parse each timing line in the hierarchical section only
    timing_lines = []
    for i in range(start_idx, end_idx):
        line = lines[i]
        if line.strip() and '[' in line and ']' in line and not line.strip().startswith('Timer report'):
            timing_lines.append(line)
    
    # Build the hierarchical structure
    root = None
    stack = []
    
    for line in timing_lines:
        parsed = parse_timing_line(line)
        if parsed:
            # Determine the level based on indentation
            level = get_indentation_level(line)
            
            # Adjust stack to current level
            while len(stack) > level:
                stack.pop()
            
            # Add to parent if exists
            if stack:
                parent = stack[-1]
                if 'children' not in parent:
                    parent['children'] = []
                parent['children'].append(parsed)
            else:
                root = parsed
            
            stack.append(parsed)
    
    return root

def parse_timing_line(line: str) -> Optional[Dict[str, Any]]:
    """
    Parse a single timing line into a dictionary.
    Expected format: [id] name time percentage count flops mflops
    """
    # Remove leading whitespace but preserve for level calculation
    stripped = line.strip()
    
    # Skip lines that don't match the expected format
    if not stripped or 'Timer report' in stripped:
        return None
    
    # Extract the bracketed ID
    id_match = re.match(r'\[(\d+)\]', stripped)
    if not id_match:
        return None
    
    timer_id = int(id_match.group(1))
    
    # Remove the ID part and parse the rest
    rest = stripped[id_match.end():].strip()
    
    # Split by whitespace
    parts = rest.split()
    
    if len(parts) < 6:
        return None
    
    # The name is everything up to the first number (time)
    # Find where the numeric data starts
    numeric_start = -1
    for i, part in enumerate(parts):
        try:
            float(part)
            numeric_start = i
            break
        except ValueError:
            continue
    
    if numeric_start == -1:
        return None
    
    name = ' '.join(parts[:numeric_start])
    numeric_parts = parts[numeric_start:]
    
    if len(numeric_parts) < 6:
        return None
    
    try:
        time_val = float(numeric_parts[0])
        percentage = float(numeric_parts[1].rstrip('%'))
        count = int(numeric_parts[2])
        flops = int(numeric_parts[3])
        mflops = float(numeric_parts[4])
        # numeric_parts[5] is "MFlops" unit
        
        return {
            'id': timer_id,
            'name': name,
            'time': time_val,
            'percentage': percentage,
            'count': count,
            'flops': flops,
            'mflops': mflops,
            'children': []
        }
    except (ValueError, IndexError):
        return None

def get_indentation_level(line: str) -> int:
    """
    Calculate the indentation level based on leading spaces.
    Each level appears to be 3 spaces in your format.
    """
    leading_spaces = len(line) - len(line.lstrip())
    return leading_spaces // 3

def print_timing_tree(node: Dict[str, Any], indent: int = 0) -> None:
    """
    Print the parsed timing tree.
    """
    if not node:
        return
    
    prefix = "  " * indent
    print(f"{prefix}[{node['id']}] {node['name']}: {node['time']:.5f}s ({node['percentage']:.1f}%) - {node['count']} calls")
    
    for child in node.get('children', []):
        print_timing_tree(child, indent + 1)

def get_timing_summary(node: Dict[str, Any]) -> Dict[str, Any]:
    """
    Get a summary of the timing data.
    """
    def collect_stats(n, stats):
        stats['total_timers'] += 1
        stats['total_time'] += n['time']
        stats['total_calls'] += n['count']
        
        if n['time'] > stats['max_time']:
            stats['max_time'] = n['time']
            stats['slowest_timer'] = n['name']
        
        for child in n.get('children', []):
            collect_stats(child, stats)
    
    stats = {
        'total_timers': 0,
        'total_time': 0,
        'total_calls': 0,
        'max_time': 0,
        'slowest_timer': ''
    }
    
    if node:
        collect_stats(node, stats)
    
    return stats

def find_timer_by_name(node: Dict[str, Any], target_name: str) -> Optional[Dict[str, Any]]:
    """
    Find the first timer with the given name in the timing tree.
    Returns the timer node or None if not found.
    """
    if not node:
        return None
    
    # Check if current node matches
    if node['name'] == target_name:
        return node
    
    # Search in children
    for child in node.get('children', []):
        result = find_timer_by_name(child, target_name)
        if result:
            return result
    
    return None

# Usage example:
if __name__ == "__main__":

    # Parse the timing report
    timing_data = parse_timing_report("""
-----------
Timer report 0 (23 timers)
--------------
[0] root 0.48496 100.0% 1 0    0.0 MFlops
   [1] BoxData::define(Box) 0.48427 99.9% 5 0    0.0 MFlops
      [2] BoxData::setVal 0.19700 40.6% 5 0    0.0 MFlops
         [3] forallInPlaceBase 0.19700 40.6% 5 0    0.0 MFlops
            [4] protoForall 0.19700 40.6% 5 0    0.0 MFlops
               [5] BoxData::makevars 0.19698 40.6% 5 0    0.0 MFlops
                  [6] indexer 0.19698 40.6% 5 0    0.0 MFlops
               [19] BoxData::var 0.00001  0.0% 5 0    0.0 MFlops
                  [21] set pointers 0.00001  0.0% 17 0    0.0 MFlops
   [7] interestingstuff 0.00065  0.1% 1 0    0.0 MFlops
      [8] forallInPlaceBase 0.00012  0.0% 10 0    0.0 MFlops
         [9] protoForall 0.00012  0.0% 10 0    0.0 MFlops
            [10] BoxData::makevars 0.00008  0.0% 10 0    0.0 MFlops
               [11] indexer 0.00008  0.0% 10 0    0.0 MFlops
            [12] BoxData::var 0.00003  0.0% 20 0    0.0 MFlops
               [16] set pointers 0.00002  0.0% 80 0    0.0 MFlops
   [13] BoxData::setVal 0.00003  0.0% 4 0    0.0 MFlops
      [14] forallInPlaceBase 0.00003  0.0% 4 0    0.0 MFlops
         [15] protoForall 0.00002  0.0% 4 0    0.0 MFlops
            [17] BoxData::makevars 0.00001  0.0% 4 0    0.0 MFlops
               [18] indexer 0.00001  0.0% 4 0    0.0 MFlops
            [20] BoxData::var 0.00001  0.0% 4 0    0.0 MFlops
               [22] set pointers 0.00000  0.0% 16 0    0.0 MFlops
    """)
    
    # # Alternatively, load from file
    # with open('TIMINGS.txt', 'r') as f:
    #     content = f.read()
    # timing_data = parse_timing_report(content)
    
    if timing_data:
        print("Parsed Timing Tree:")
        print_timing_tree(timing_data)
        
        print("\n" + "="*50 + "\n")
        print("Summary:")
        summary = get_timing_summary(timing_data)
        print(f"Total timers: {summary['total_timers']}")
        print(f"Root execution time: {timing_data['time']:.5f}s")
        print(f"Total calls: {summary['total_calls']}")
        print(f"Slowest operation: {summary['slowest_timer']} ({summary['max_time']:.5f}s)")
    else:
        print("Failed to parse timing data")

    # Get information for a specific node
    interestingstuff_timer = find_timer_by_name(timing_data, 'interestingstuff')

    if interestingstuff_timer:
        print("\n" + "="*50 + "\n")
        print("Get info for a specific node (e.g., 'root'):")
        print(f"Found {interestingstuff_timer['name']}:")
        print(f"  Time: {interestingstuff_timer['time']:.6f}s")
        print(f"  Percentage: {interestingstuff_timer['percentage']:.1f}%")
        print(f"  Calls: {interestingstuff_timer['count']}")
        print(f"  Children: {len(interestingstuff_timer['children'])}")
    else:
        print("Timer 'interestingstuff' not found")