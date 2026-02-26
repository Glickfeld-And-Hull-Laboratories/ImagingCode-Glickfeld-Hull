#!/usr/bin/env python3
"""
Combine step3_poolSingleDay PDF outputs into a single consolidated report.
Creates a summary page with experimental details and combines all plots.

Usage: python combine_step3_report.py /path/to/step3/output/folder
"""

import os
import sys
import re
from pathlib import Path
import argparse
from datetime import datetime

try:
    from reportlab.pdfgen import canvas
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
    from reportlab.lib.units import inch
    from PyPDF2 import PdfMerger
except ImportError:
    print("Required packages missing. Install with:")
    print("pip install reportlab PyPDF2")
    sys.exit(1)


def parse_folder_name(folder_path):
    """Extract experimental details from step3 output folder structure."""
    folder_name = os.path.basename(folder_path)
    parent_name = os.path.basename(os.path.dirname(folder_path))
    
    # Parse session list from parent folder (e.g., "sess1_2_3")
    sess_match = re.match(r'sess(.+)', parent_name)
    if sess_match:
        sessions = sess_match.group(1).replace('_', ', ')
    else:
        sessions = "Unknown"
    
    # Parse date and filter info from folder name
    date_match = re.search(r'(\d{4}-\d{2}-\d{2})', folder_name)
    date_str = date_match.group(1) if date_match else "Unknown date"
    
    # Parse size filters
    red_filter = None
    green_filter = None
    
    red_match = re.search(r'_red([\d.]+)deg', folder_name)
    if red_match:
        red_filter = float(red_match.group(1))
        
    green_match = re.search(r'_green([\d.]+)deg', folder_name)
    if green_match:
        green_filter = float(green_match.group(1))
    
    return {
        'sessions': sessions,
        'date': date_str,
        'red_filter': red_filter,
        'green_filter': green_filter
    }


def create_summary_page(output_path, folder_info, pdf_files):
    """Create a summary page with experimental details."""
    doc = SimpleDocTemplate(output_path, pagesize=letter)
    styles = getSampleStyleSheet()
    story = []
    
    # Title
    title = Paragraph("Pooled Single-Day Analysis Report", styles['Title'])
    story.append(title)
    story.append(Spacer(1, 0.3*inch))
    
    # Analysis details
    story.append(Paragraph("Analysis Details", styles['Heading2']))
    story.append(Spacer(1, 0.1*inch))
    
    details = [
        f"<b>Sessions analyzed:</b> {folder_info['sessions']}",
        f"<b>Analysis date:</b> {folder_info['date']}",
        f"<b>Generated on:</b> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        ""
    ]
    
    # Size filter information
    if folder_info['red_filter'] and folder_info['green_filter']:
        details.append(f"<b>Size filters applied:</b>")
        details.append(f"• RED (HTP+) cells: responsive at {folder_info['red_filter']}°")
        details.append(f"• GREEN (HTP-) cells: responsive at {folder_info['green_filter']}°")
    elif folder_info['red_filter']:
        details.append(f"<b>Size filter applied:</b> RED (HTP+) cells responsive at {folder_info['red_filter']}° (GREEN cells unfiltered)")
    elif folder_info['green_filter']:
        details.append(f"<b>Size filter applied:</b> GREEN (HTP-) cells responsive at {folder_info['green_filter']}° (RED cells unfiltered)")
    else:
        details.append("<b>Size filters:</b> None (all responsive cells included)")
    
    details.extend([
        "",
        "<b>Contents:</b>",
        "• Size tuning curves - stationary conditions",
        "• Neural timecourses - stationary conditions", 
        "• Size tuning curves - stationary vs running (condition-matched cells)",
        "• Neural timecourses - stationary vs running (condition-matched cells)",
        "• Bonus: Size tuning with contrast overlays",
        "",
        f"<b>Total figures:</b> {len(pdf_files)} plots"
    ])
    
    for detail in details:
        if detail:
            story.append(Paragraph(detail, styles['Normal']))
        else:
            story.append(Spacer(1, 0.1*inch))
    
    doc.build(story)


def get_plot_order():
    """Define the logical order for combining plots."""
    return [
        ('sizeResponse_stat.pdf', 'Size Tuning - Stationary'),
        ('timecourse_stat_', 'Neural Timecourses - Stationary'),
        ('sizeResponse_byCondition.pdf', 'Size Tuning - Stationary vs Running'),
        ('timecourse_byCondition_', 'Neural Timecourses - Stationary vs Running'),
        ('sizeResponse_stat_contrastOverlay.pdf', 'Size Tuning - Contrast Overlay (Stationary)'),
        ('sizeResponse_contrastOverlay_byCondition.pdf', 'Size Tuning - Contrast Overlay (Stat vs Run)')
    ]


def find_matching_files(pdf_files, pattern):
    """Find all files matching a pattern (for timecourse files with multiple sizes)."""
    if pattern.endswith('_'):
        # Pattern like 'timecourse_stat_' - find all matching files
        matching = [f for f in pdf_files if f.startswith(pattern)]
        return sorted(matching)  # Sort to ensure consistent order
    else:
        # Exact filename match
        return [pattern] if pattern in pdf_files else []


def combine_pdfs(folder_path, output_name=None):
    """Combine all PDFs in the folder into a single report."""
    folder_path = Path(folder_path)
    
    if not folder_path.exists():
        print(f"Error: Folder {folder_path} does not exist")
        return False
    
    # Find all PDF files
    pdf_files = [f.name for f in folder_path.glob('*.pdf')]
    
    if not pdf_files:
        print(f"Error: No PDF files found in {folder_path}")
        return False
    
    # Parse folder information
    folder_info = parse_folder_name(folder_path)
    
    # Create output filename
    if output_name is None:
        base_name = f"pooled_analysis_sess{folder_info['sessions'].replace(', ', '_')}"
        if folder_info['red_filter'] and folder_info['green_filter']:
            base_name += f"_red{folder_info['red_filter']}_green{folder_info['green_filter']}"
        elif folder_info['red_filter']:
            base_name += f"_red{folder_info['red_filter']}"
        elif folder_info['green_filter']:
            base_name += f"_green{folder_info['green_filter']}"
        output_name = f"{base_name}.pdf"
    
    output_path = folder_path / output_name
    temp_summary = folder_path / "temp_summary.pdf"
    
    try:
        # Create summary page
        print("Creating summary page...")
        create_summary_page(str(temp_summary), folder_info, pdf_files)
        
        # Initialize PDF merger
        merger = PdfMerger()
        
        # Add summary page first
        merger.append(str(temp_summary))
        
        # Add plots in logical order
        plot_order = get_plot_order()
        added_files = set()
        
        print("Adding plots in logical order...")
        for pattern, description in plot_order:
            matching_files = find_matching_files(pdf_files, pattern)
            for pdf_file in matching_files:
                if pdf_file not in added_files:
                    pdf_path = folder_path / pdf_file
                    print(f"  Adding: {pdf_file}")
                    merger.append(str(pdf_path))
                    added_files.add(pdf_file)
        
        # Add any remaining PDFs that weren't caught by the patterns
        remaining_files = set(pdf_files) - added_files
        if remaining_files:
            print("Adding remaining files:")
            for pdf_file in sorted(remaining_files):
                pdf_path = folder_path / pdf_file
                print(f"  Adding: {pdf_file}")
                merger.append(str(pdf_path))
        
        # Write combined PDF
        print(f"Writing combined PDF: {output_path}")
        merger.write(str(output_path))
        merger.close()
        
        # Clean up temporary summary file
        temp_summary.unlink()
        
        print(f"✓ Successfully created: {output_path}")
        print(f"  Combined {len(pdf_files)} plots with summary page")
        
        return True
        
    except Exception as e:
        print(f"Error combining PDFs: {e}")
        # Clean up temp file if it exists
        if temp_summary.exists():
            temp_summary.unlink()
        return False


def main():
    parser = argparse.ArgumentParser(description='Combine step3_poolSingleDay PDF outputs')
    parser.add_argument('folder', help='Path to step3 output folder containing PDFs')
    parser.add_argument('-o', '--output', help='Output filename (optional)')
    
    args = parser.parse_args()
    
    success = combine_pdfs(args.folder, args.output)
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
