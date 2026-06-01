# -*- coding: utf-8 -*-
"""
created by GROUP 78 MEMBERS 
Continuous beam solver using  three-moment equations.

Sign convention:
    - Downward loads are entered as positive values.
    - Sagging bending moments are positive.
    - Hogging support moments normally appear as negative values.
"""

import json
import tkinter as tk
from tkinter import Canvas, messagebox, simpledialog, ttk


app = None
menuBar = None

# GLOBAL VARIABLES

Support = []
ListofSpans = []

EndSupportMoments = {"Left Moment": 0.0,"Right Moment": 0.0,}

Overhangs = {"Left": {"enabled": False, "Length": 0.0, "Load Type": "None"},"Right": {"enabled": False, "Length": 0.0, "Load Type": "None"},}

SolvedMoments = []
SolvedReactions = []
AnalysisResult = {}
EquationRows = []
LoadInputFrame = None


# SMALL UTILITIES

def safe_float(value, default=0.0):
    try:
        if value is None or value == "":
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def clamp(value, low, high):
    return max(low, min(high, value))


def fmt(value, digits=4):
    value = safe_float(value)
    if abs(value) < 10 ** (-(digits + 1)):
        value = 0.0
    text = f"{value:.{digits}g}"
    return "0" if text == "-0" else text


def show_error(title, message):
    if app is not None:
        messagebox.showerror(title, message)
    print(f"{title}: {message}")


def show_info(title, message):
    if app is not None:
        messagebox.showinfo(title, message)
    print(f"{title}: {message}")


def write_json(filename, data):
    with open(filename, "w", encoding="utf-8") as file:
        json.dump(data, file, indent=4)


def solve_linear_system(matrix, rhs):
    """Solve a dense linear system with partial-pivot Gaussian elimination."""

    n = len(rhs)
    aug = [[float(matrix[i][j]) for j in range(n)] + [float(rhs[i])]
        for i in range(n)]

    for col in range(n):
        pivot = max(range(col, n), key=lambda row: abs(aug[row][col]))
        if abs(aug[pivot][col]) < 1e-12:
            raise ValueError("The three-moment matrix is singular. Check supports, spans, and EI values.")

        if pivot != col:
            aug[col], aug[pivot] = aug[pivot], aug[col]

        pivot_value = aug[col][col]
        for j in range(col, n + 1):
            aug[col][j] /= pivot_value

        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            if abs(factor) < 1e-14:
                continue
            for j in range(col, n + 1):
                aug[row][j] -= factor * aug[col][j]

    return [aug[row][n] for row in range(n)]


def support_type(index):
    if 0 <= index < len(Support):
        return Support[index].get("Support Type", "")
    return ""


def support_is_fixed(index):
    return support_type(index) == "Fixed Support"


def span_length(span):
    return safe_float(span.get("Span Length"), 0.0)


def span_ei(span):
    ei = safe_float(span.get("EI"), 1.0)
    return ei if ei > 0 else 1.0


def validate_model():
    errors = []

    if not ListofSpans:
        errors.append("Enter at least one span.")

    if len(Support) != len(ListofSpans) + 1:
        errors.append("Number of supports must be one more than the number of spans.")

    for index, span in enumerate(ListofSpans):
        if not span:
            errors.append(f"Span {index + 1} has not been saved.")
            continue

        if span_length(span) <= 0:
            errors.append(f"Span {index + 1} length must be greater than zero.")

        if safe_float(span.get("EI"), 1.0) <= 0:
            errors.append(f"Span {index + 1} EI must be greater than zero.")

    for index, support in enumerate(Support):
        if not support:
            errors.append(f"Support {index + 1} has not been saved.")
            continue

        if support.get("Support Type") not in ["Fixed Support","Hinge Support","Roller Support",]:
            errors.append(f"Support {index + 1} type is missing.")

    if errors:
        raise ValueError("\n".join(errors))


#  LOAD NORMALIZATION 

def normalized_span_loads(span):
    """Return loads on a main span in a common format."""

    L = span_length(span)
    load_type = span.get("Load Type", "No Load")
    loads = []

    if load_type == "Point Load":
        P = safe_float(span.get("Point Load Value"), 0.0)
        a = clamp(safe_float(span.get("Point Load Position"), 0.0), 0.0, L)
        if P != 0:
            loads.append({"kind": "point","P": P,"x": a,"W": P,"centroid": a,})

    elif load_type == "UDL":
        w = safe_float(span.get("UDL Intensity"), 0.0)
        start = clamp(safe_float(span.get("UDL Start"), 0.0), 0.0, L)
        length = safe_float(span.get("UDL Length"), L - start)
        end = clamp(start + max(length, 0.0), start, L)
        length = end - start

        if w != 0 and length > 0:
            loads.append({"kind": "udl","w": w,"start": start,"end": end,"length": length,"W": w * length,"centroid": start + length / 2.0,})

    return loads


def load_resultants(loads, L):
    total_load = sum(load["W"] for load in loads)
    moment_about_left = sum(load["W"] * load["centroid"] for load in loads)
    moment_about_right = sum(load["W"] * (L - load["centroid"]) for load in loads)
    return total_load, moment_about_left, moment_about_right


def simple_span_load_stats(span):
    """Exact area and first moments of the simply-supported BMD caused by span loads.
       These are the A*xbar terms used in the three-moment equation:
        Ax_left  = integral( x * M0(x) dx )
        Ax_right = integral( (L - x) * M0(x) dx )
    """

    L = span_length(span)
    loads = normalized_span_loads(span)
    total_load, moment_about_left, moment_about_right = load_resultants(loads, L)

    if L <= 0:
        return {"A": 0.0,"Ax_left": 0.0,"Ax_right": 0.0,"xbar_left": 0.0,"xbar_right": 0.0,"simple_R_left": 0.0,"simple_R_right": 0.0,"total_load": total_load,}

    R_left = moment_about_right / L
    R_right = total_load - R_left

    area = R_left * L ** 2 / 2.0
    first_left = R_left * L ** 3 / 3.0

    for load in loads:
        if load["kind"] == "point":
            P = load["P"]
            a = load["x"]
            D = L - a
            area -= P * D ** 2 / 2.0
            first_left -= P * (D ** 3 / 3.0 + a * D ** 2 / 2.0)

        elif load["kind"] == "udl":
            w = load["w"]
            start = load["start"]
            end = load["end"]
            Ds = L - start
            De = L - end

            area -= w * (Ds ** 3 - De ** 3) / 6.0

            Js = Ds ** 4 / 4.0 + start * Ds ** 3 / 3.0
            Je = De ** 4 / 4.0 + end * De ** 3 / 3.0
            first_left -= w * (Js - Je) / 2.0

    first_right = L * area - first_left

    if abs(area) > 1e-12:
        xbar_left = first_left / area
        xbar_right = first_right / area
    else:
        xbar_left = 0.0
        xbar_right = 0.0

    return {"A": area,"Ax_left": first_left,"Ax_right": first_right,"xbar_left": xbar_left,"xbar_right": xbar_right,"simple_R_left": R_left,
            "simple_R_right": R_right,"total_load": total_load,"moment_about_left": moment_about_left,"moment_about_right": moment_about_right,}


#  OVERHANGS

def normalized_overhang_loads(side):
    """Loads are measured from the support toward the free end."""

    data = Overhangs.get(side, {})
    if not data.get("enabled"):
        return []

    L = safe_float(data.get("Length"), 0.0)
    if L <= 0:
        return []

    load_type = data.get("Load Type", "None")
    loads = []

    if load_type == "Point Load":
        P = safe_float(data.get("Point Load Value"), 0.0)
        d = clamp(safe_float(data.get("Point Load Distance"), L), 0.0, L)
        if P != 0:
            loads.append({"kind": "point","P": P,"d": d,"W": P,"centroid_from_support": d,})

    elif load_type == "UDL":
        w = safe_float(data.get("UDL Intensity"), 0.0)
        start = clamp(safe_float(data.get("UDL Start"), 0.0), 0.0, L)
        length = safe_float(data.get("UDL Length"), L - start)
        end = clamp(start + max(length, 0.0), start, L)
        length = end - start

        if w != 0 and length > 0:
            loads.append({"kind": "udl","w": w,"start": start,"end": end,"length": length,"W": w * length,"centroid_from_support": start + length / 2.0,})

    return loads


def overhang_length(side):
    data = Overhangs.get(side, {})
    return safe_float(data.get("Length"), 0.0) if data.get("enabled") else 0.0


def overhang_resultant(side):
    loads = normalized_overhang_loads(side)
    total = sum(load["W"] for load in loads)
    moment_at_support = -sum(
        load["W"] * load["centroid_from_support"] for load in loads)
    return total, moment_at_support


def boundary_moment_from_end_condition(side):
    key = "Left Moment" if side == "Left" else "Right Moment"
    direct_moment = safe_float(EndSupportMoments.get(key), 0.0)
    _, overhang_moment = overhang_resultant(side)
    return direct_moment + overhang_moment


def ask_yes_no(title, prompt, default="no"):
    answer = simpledialog.askstring(title, prompt, initialvalue=default)
    if answer is None:
        return default.lower().startswith("y")
    return answer.strip().lower().startswith("y")


def get_end_support_moments():
    """Collect direct end moments and real overhang load data."""

    global EndSupportMoments, Overhangs

    for side in ["Left", "Right"]:
        key = f"{side} Moment"
        current_direct = safe_float(EndSupportMoments.get(key), 0.0)

        direct = simpledialog.askfloat(
            f"{side} End Moment",
            (f"Enter direct {side.lower()} end moment in kNm.\n"
                "Use positive for sagging and negative for hogging."),initialvalue=current_direct,)

        if direct is not None:
            EndSupportMoments[key] = direct

        has_overhang = ask_yes_no(f"{side} Overhang", f"Is there a real {side.lower()} overhang span? (yes/no)",
            "yes" if Overhangs.get(side, {}).get("enabled") else "no",)

        if not has_overhang:
            Overhangs[side] = {"enabled": False, "Length": 0.0, "Load Type": "None"}
            continue

        length = simpledialog.askfloat(f"{side} Overhang Length","Enter overhang length in m.",initialvalue=safe_float(Overhangs.get(side, {}).get("Length"), 1.0),)
        if length is None or length <= 0:
            Overhangs[side] = {"enabled": False, "Length": 0.0, "Load Type": "None"}
            continue

        load_type = simpledialog.askstring(f"{side} Overhang Load","Enter load type: point, udl, or none.",initialvalue=Overhangs.get(side, {}).get("Load Type", "point"),)

        load_type = (load_type or "none").strip().lower()

        if load_type == "point":
            force = simpledialog.askfloat(f"{side} Overhang Point Load","Enter point load value in kN.",initialvalue=safe_float(Overhangs.get(side, {}).get("Point Load Value"), 0.0),)
            distance = simpledialog.askfloat(f"{side} Overhang Point Load","Enter load distance from the support toward the free end in m.",initialvalue=safe_float(Overhangs.get(side, {}).get("Point Load Distance"), length),)

            Overhangs[side] = {"enabled": True,"Length": length,"Load Type": "Point Load","Point Load Value": force or 0.0,"Point Load Distance": distance if distance is not None else length,}

        elif load_type == "udl":
            intensity = simpledialog.askfloat(f"{side} Overhang UDL","Enter UDL intensity in kN/m.",initialvalue=safe_float(Overhangs.get(side, {}).get("UDL Intensity"), 0.0),)
            start = simpledialog.askfloat(f"{side} Overhang UDL","Enter UDL start distance from the support in m.", initialvalue=safe_float(Overhangs.get(side, {}).get("UDL Start"), 0.0),)
            udl_length = simpledialog.askfloat(f"{side} Overhang UDL","Enter loaded length in m.",initialvalue=safe_float(Overhangs.get(side, {}).get("UDL Length"), length ),
            )

            Overhangs[side] = {"enabled": True,"Length": length,"Load Type": "UDL","UDL Intensity": intensity or 0.0,"UDL Start": start or 0.0,"UDL Length": udl_length if udl_length is not None else length, }

        else:
            Overhangs[side] = {"enabled": True,"Length": length,"Load Type": "None",}

    write_json("EndConditions.json",{"EndSupportMoments": EndSupportMoments, "Overhangs": Overhangs},)


# THREE-MOMENT SOLVER

def make_numeric_equation(coeffs, rhs):
    terms = []
    for index, coeff in enumerate(coeffs):
        if abs(coeff) > 1e-12:
            terms.append(f"({fmt(coeff, 5)})M{index + 1}")
    if not terms:
        terms.append("0")
    return " + ".join(terms).replace("+ (-", "- (") + f" = {fmt(rhs, 6)}"


def analyze_beam():
    global SolvedMoments, SolvedReactions, AnalysisResult, EquationRows

    validate_model()

    spans = ListofSpans
    span_count = len(spans)
    support_count = span_count + 1
    stats = [simple_span_load_stats(span) for span in spans]

    rows = []

    # Left boundary equation.
    left_row = [0.0] * support_count
    if support_is_fixed(0):
        L = span_length(spans[0])
        EI = span_ei(spans[0])
        left_row[0] = 2.0 * L / EI
        left_row[1] = L / EI
        rhs = -6.0 * stats[0]["Ax_right"] / (L * EI)
        symbol = (f"Fixed left: 2({fmt(L)}/{fmt(EI)})M1 + "f"({fmt(L)}/{fmt(EI)})M2 = "f"-6(A1*xbar_right)/({fmt(L)}*{fmt(EI)}) = {fmt(rhs, 6)}")
    else:
        left_row[0] = 1.0
        rhs = boundary_moment_from_end_condition("Left")
        symbol = f"Simple/overhang left: M1 = {fmt(rhs, 6)}"

    rows.append({"name": "Left boundary","coeffs": left_row,"rhs": rhs,"symbol": symbol,})

    # Interior three-moment equations.
    for i in range(span_count - 1):
        span1 = spans[i]
        span2 = spans[i + 1]
        L1 = span_length(span1)
        L2 = span_length(span2)
        EI1 = span_ei(span1)
        EI2 = span_ei(span2)
        S1 = stats[i]
        S2 = stats[i + 1]

        row = [0.0] * support_count
        row[i] = L1 / EI1
        row[i + 1] = 2.0 * (L1 / EI1 + L2 / EI2)
        row[i + 2] = L2 / EI2

        rhs = -6.0 * (S1["Ax_left"] / (L1 * EI1)+ S2["Ax_right"] / (L2 * EI2))

        symbol = (f"Span {i + 1}-{i + 2}: " f"({fmt(L1)}/{fmt(EI1)})M{i + 1} + "f"2[({fmt(L1)}/{fmt(EI1)}) + ({fmt(L2)}/{fmt(EI2)})]M{i + 2} + "
                  f"({fmt(L2)}/{fmt(EI2)})M{i + 3} = "f"-6[(A{i + 1}*xbar_left)/({fmt(L1)}*{fmt(EI1)}) + "f"(A{i + 2}*xbar_right)/({fmt(L2)}*{fmt(EI2)})] = {fmt(rhs, 6)}")

        rows.append({"name": f"Interior support {i + 2}","coeffs": row,"rhs": rhs,"symbol": symbol,})

    # Right boundary equation.
    right_row = [0.0] * support_count
    if support_is_fixed(support_count - 1):
        L = span_length(spans[-1])
        EI = span_ei(spans[-1])
        right_row[-2] = L / EI
        right_row[-1] = 2.0 * L / EI
        rhs = -6.0 * stats[-1]["Ax_left"] / (L * EI)
        symbol = (f"Fixed right: ({fmt(L)}/{fmt(EI)})M{support_count - 1} + "f"2({fmt(L)}/{fmt(EI)})M{support_count} = "f"-6(A{span_count}*xbar_left)/({fmt(L)}*{fmt(EI)}) = {fmt(rhs, 6)}")
    else:
        right_row[-1] = 1.0
        rhs = boundary_moment_from_end_condition("Right")
        symbol = f"Simple/overhang right: M{support_count} = {fmt(rhs, 6)}"

    rows.append({"name": "Right boundary","coeffs": right_row,"rhs": rhs,"symbol": symbol,})

    A = [row["coeffs"][:] for row in rows]
    B = [row["rhs"] for row in rows]
    moments = solve_linear_system(A, B)

    span_forces = []
    reactions = [0.0] * support_count

    for index, span in enumerate(spans):
        L = span_length(span)
        loads = normalized_span_loads(span)
        total_load, _, moment_about_right = load_resultants(loads, L)
        M_left = moments[index]
        M_right = moments[index + 1]

        R_left = (M_right - M_left + moment_about_right) / L
        R_right = total_load - R_left

        reactions[index] += R_left
        reactions[index + 1] += R_right

        span_forces.append({"R_left": R_left,"R_right": R_right,"M_left": M_left,"M_right": M_right,"loads": loads,})

    left_overhang_load, _ = overhang_resultant("Left")
    right_overhang_load, _ = overhang_resultant("Right")
    if reactions:
        reactions[0] += left_overhang_load
        reactions[-1] += right_overhang_load

    SolvedMoments = [round(float(value), 6) for value in moments]
    SolvedReactions = [round(float(value), 6) for value in reactions]
    EquationRows = rows
    AnalysisResult = {
        "moments": moments,
        "reactions": reactions,
        "span_forces": span_forces,
        "span_stats": stats,
        "equations": rows,
        "matrix_A": A,
        "matrix_B": B,
    }

    return AnalysisResult


def three_moment_equation():
    try:
        if app is not None and messagebox.askyesno(
            "End Conditions",
            "Edit direct end moments and overhangs before solving?",
        ):
            get_end_support_moments()

        result = analyze_beam()

    except ValueError as exc:
        show_error("Beam Solver", str(exc))
        return

    print("\n========== THREE-MOMENT EQUATIONS ==========\n")
    for row in result["equations"]:
        print(row["symbol"])
        print(make_numeric_equation(row["coeffs"], row["rhs"]))

    print("\n========== SUPPORT MOMENTS ==========")
    for i, moment in enumerate(SolvedMoments):
        print(f"M{i + 1} = {moment} kNm")

    print("\n========== SUPPORT REACTIONS ==========")
    for i, reaction in enumerate(SolvedReactions):
        print(f"R{i + 1} = {reaction} kN")

    show_equations()


# SHEAR AND MOMENT FUNCTIONS 

def span_shear_moment(span, x, M_left, R_left):
    L = span_length(span)
    x = clamp(x, 0.0, L)
    V = R_left
    M = M_left + R_left * x

    for load in normalized_span_loads(span):
        if load["kind"] == "point":
            a = load["x"]
            if x >= a:
                V -= load["P"]
                M -= load["P"] * (x - a)

        elif load["kind"] == "udl":
            start = load["start"]
            end = load["end"]

            if x > start:
                loaded_length = min(x, end) - start
                loaded_length = max(loaded_length, 0.0)
                V -= load["w"] * loaded_length

                M -= load["w"] * loaded_length * (x - (start + loaded_length / 2.0))

    return V, M


def overhang_samples(side, count=80):
    L = overhang_length(side)
    if L <= 0:
        return []

    loads = normalized_overhang_loads(side)
    samples = []

    if side == "Left":
        # Local x is measured from the free end toward the first support.
        for k in range(count + 1):
            x = L * k / count
            V = 0.0
            M = 0.0

            for load in loads:
                if load["kind"] == "point":
                    load_x = L - load["d"]
                    if x >= load_x:
                        V -= load["P"]
                        M -= load["P"] * (x - load_x)

                elif load["kind"] == "udl":
                    load_start = L - load["end"]
                    load_end = L - load["start"]
                    if x > load_start:
                        loaded_length = min(x, load_end) - load_start
                        loaded_length = max(loaded_length, 0.0)
                        V -= load["w"] * loaded_length
                        M -= load["w"] * loaded_length * (x - (load_start + loaded_length / 2.0))

            samples.append((x, V, M))

    else:
        # Local x is measured from the last support toward the free end.
        total_load, support_moment = overhang_resultant("Right")
        for k in range(count + 1):
            x = L * k / count
            V = total_load
            M = support_moment + total_load * x

            for load in loads:
                if load["kind"] == "point":
                    if x >= load["d"]:
                        V -= load["P"]
                        M -= load["P"] * (x - load["d"])

                elif load["kind"] == "udl":
                    if x > load["start"]:
                        loaded_length = min(x, load["end"]) - load["start"]
                        loaded_length = max(loaded_length, 0.0)
                        V -= load["w"] * loaded_length
                        M -= load["w"] * loaded_length * (x - (load["start"] + loaded_length / 2.0))

            samples.append((x, V, M))

    return samples


def load_breakpoints(span):
    L = span_length(span)
    points = {0.0, L}

    for load in normalized_span_loads(span):
        if load["kind"] == "point":
            points.add(load["x"])
        elif load["kind"] == "udl":
            points.add(load["start"])
            points.add(load["end"])

    extra = set()
    eps = max(L * 1e-6, 1e-7)
    for point in points:
        extra.add(clamp(point - eps, 0.0, L))
        extra.add(point)
        extra.add(clamp(point + eps, 0.0, L))

    for k in range(1, 60):
        extra.add(L * k / 60.0)

    return sorted(extra)


def beam_geometry():
    left = overhang_length("Left")
    right = overhang_length("Right")
    support_positions = [left]
    x = left

    for span in ListofSpans:
        x += span_length(span)
        support_positions.append(x)

    return {"left_overhang": left,"right_overhang": right,"support_positions": support_positions,"total_length": x + right,}


def diagram_points(kind):
    if not AnalysisResult:
        analyze_beam()

    geometry = beam_geometry()
    points = []

    for x, V, M in overhang_samples("Left"):
        value = V if kind == "SFD" else M
        points.append((x, value))

    left_offset = geometry["left_overhang"]
    x_start = left_offset

    for index, span in enumerate(ListofSpans):
        forces = AnalysisResult["span_forces"][index]
        for local_x in load_breakpoints(span):
            V, M = span_shear_moment(span,local_x,forces["M_left"],forces["R_left"],)
            value = V if kind == "SFD" else M
            points.append((x_start + local_x, value))

        x_start += span_length(span)

    right_samples = overhang_samples("Right")
    if right_samples:
        right_start = geometry["support_positions"][-1]
        for x, V, M in right_samples:
            value = V if kind == "SFD" else M
            points.append((right_start + x, value))

    return sorted(points, key=lambda item: item[0])


# TKINTER MODEL INPUT 
def newModel():
    global Support, ListofSpans, LoadInputFrame

    if app is None:
        return

    Support = []
    ListofSpans = []
    LoadInputFrame = None

    try:
        Beamprpts.destroy()
    except NameError:
        pass
    except tk.TclError:
        pass

    beam_frame = tk.LabelFrame(app, text="Beam Properties")
    beam_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nw")
    globals()["Beamprpts"] = beam_frame

    label1 = tk.Label(beam_frame, text="Number of supports")
    label1.grid(row=0, column=0, sticky="w")

    NumberofSupports = tk.Entry(beam_frame)
    NumberofSupports.grid(row=0, column=1)

    def Support_Proceed():
        global Support, comboSupport

        try:
            value = int(NumberofSupports.get())
        except ValueError:
            show_error("Supports", "Enter a valid number of supports.")
            return

        if value < 2:
            show_error("Supports", "A beam needs at least two supports.")
            return

        Support = [{} for _ in range(value)]
        listOfSupportNames = [f"Support {i + 1}" for i in range(value)]

        comboSupport = ttk.Combobox(beam_frame, values=listOfSupportNames)
        comboSupport.grid(row=0, column=3)
        comboSupport.bind("<<ComboboxSelected>>", support_selected)

    def support_selected(event):
        try:
            inputFrame.destroy()
        except NameError:
            pass
        except tk.TclError:
            pass

        frame = tk.LabelFrame(app, text=f"{comboSupport.get()} Information")
        frame.grid(row=2, column=1, padx=10, pady=10, sticky="nw")
        globals()["inputFrame"] = frame

        Supp_Type = ttk.Label(frame, text="Select Support Type")
        Supp_Type.grid(row=0, column=0)

        Supp_List = ttk.Combobox(frame,values=["Fixed Support", "Hinge Support", "Roller Support"],)
        Supp_List.grid(row=0, column=1)

        support_index = int(comboSupport.get().split(" ")[1]) - 1
        existing = Support[support_index] if support_index < len(Support) else {}
        if existing.get("Support Type"):
            Supp_List.set(existing["Support Type"])

        def save_support():
            Support[support_index] = {"Support Name": comboSupport.get(),"Support Type": Supp_List.get(),}
            write_json("SupportData.json", Support)
            print(Support)

        SaveButton = ttk.Button(frame, text="Save", command=save_support)
        SaveButton.grid(row=2, column=0)

    SupportProceed = ttk.Button(beam_frame,text="Proceed",command=Support_Proceed,)
    SupportProceed.grid(row=0, column=2)

    label2 = tk.Label(beam_frame, text="Number of Spans")
    label2.grid(row=1, column=0, sticky="w")

    NumberofMembers = tk.Entry(beam_frame)
    NumberofMembers.grid(row=1, column=1)

    def SpanProceed():
        global ListofSpans, comboSpanName

        try:
            value = int(NumberofMembers.get())
        except ValueError:
            show_error("Spans", "Enter a valid number of spans.")
            return

        if value < 1:
            show_error("Spans", "Enter at least one span.")
            return

        ListofSpans = [{} for _ in range(value)]
        SpanName = [f"Span {i + 1}" for i in range(value)]

        comboSpanName = ttk.Combobox(beam_frame, values=SpanName)
        comboSpanName.grid(row=1, column=3)
        comboSpanName.bind("<<ComboboxSelected>>", Span_Selected)

    def Span_Selected(event):
        global LoadInputFrame

        try:
            SimpleBeam.destroy()
        except NameError:
            pass
        except tk.TclError:
            pass

        if LoadInputFrame is not None:
            try:
                LoadInputFrame.destroy()
            except tk.TclError:
                pass
            LoadInputFrame = None

        frame = tk.LabelFrame(app, text="Span Information")
        frame.grid(row=2, column=0, padx=10, pady=10, sticky="nw")
        globals()["SimpleBeam"] = frame

        span_index = int(comboSpanName.get().split(" ")[1]) - 1
        existing = ListofSpans[span_index] if span_index < len(ListofSpans) else {}

        LT = ttk.Label(frame, text="Loading Type")
        LT.grid(row=0, column=0, sticky="w")

        comboLoad = ttk.Combobox(frame, values=["No Load", "Point Load", "UDL"])
        comboLoad.grid(row=0, column=1)
        comboLoad.set(existing.get("Load Type", "No Load"))

        SL = ttk.Label(frame, text="Span Length (m)")
        SL.grid(row=1, column=0, sticky="w")

        sl = ttk.Entry(frame)
        sl.grid(row=1, column=1)
        sl.insert(0, existing.get("Span Length", ""))

        EILabel = ttk.Label(frame, text="EI for this span")
        EILabel.grid(row=2, column=0, sticky="w")

        ei_entry = ttk.Entry(frame)
        ei_entry.grid(row=2, column=1)
        ei_entry.insert(0, existing.get("EI", "1"))

        load_entries = {}
        load_frame = {"widget": None}

        def load_selected(event=None):
            global LoadInputFrame

            if load_frame["widget"] is not None:
                load_frame["widget"].destroy()
                LoadInputFrame = None

            load_entries.clear()
            loadType = comboLoad.get()

            if loadType == "No Load":
                return

            lf = ttk.LabelFrame(app, text=f"{loadType} Properties")
            lf.grid(row=3, column=0, padx=10, pady=10, sticky="nw")
            load_frame["widget"] = lf
            LoadInputFrame = lf

            if loadType == "Point Load":
                lblPosition = ttk.Label(lf, text="Position from left (m)")
                lblPosition.grid(row=0, column=0, sticky="w")
                entPosition = ttk.Entry(lf)
                entPosition.grid(row=0, column=1)
                entPosition.insert(0, existing.get("Point Load Position", ""))

                lblValue = ttk.Label(lf, text="Point Load Value (kN)")
                lblValue.grid(row=1, column=0, sticky="w")
                entValue = ttk.Entry(lf)
                entValue.grid(row=1, column=1)
                entValue.insert(0, existing.get("Point Load Value", ""))

                load_entries["position"] = entPosition
                load_entries["value"] = entValue

            elif loadType == "UDL":
                lblIntensity = ttk.Label(lf, text="UDL Intensity (kN/m)")
                lblIntensity.grid(row=0, column=0, sticky="w")
                entIntensity = ttk.Entry(lf)
                entIntensity.grid(row=0, column=1)
                entIntensity.insert(0, existing.get("UDL Intensity", ""))

                lblStart = ttk.Label(lf, text="UDL Start from left (m)")
                lblStart.grid(row=1, column=0, sticky="w")
                entStart = ttk.Entry(lf)
                entStart.grid(row=1, column=1)
                entStart.insert(0, existing.get("UDL Start", "0"))

                lblLength = ttk.Label(lf, text="Length of UDL (m)")
                lblLength.grid(row=2, column=0, sticky="w")
                entLength = ttk.Entry(lf)
                entLength.grid(row=2, column=1)
                entLength.insert(0, existing.get("UDL Length", ""))

                load_entries["intensity"] = entIntensity
                load_entries["start"] = entStart
                load_entries["length"] = entLength

        comboLoad.bind("<<ComboboxSelected>>", load_selected)
        load_selected()

        def save_span():
            loadType = comboLoad.get()
            span_data = {"Span Name": comboSpanName.get(),"Span Length": sl.get(),"EI": ei_entry.get() or "1", "Load Type": loadType,}

            if loadType == "Point Load":
                span_data["Point Load Position"] = load_entries["position"].get()
                span_data["Point Load Value"] = load_entries["value"].get()

            elif loadType == "UDL":
                span_data["UDL Intensity"] = load_entries["intensity"].get()
                span_data["UDL Start"] = load_entries["start"].get()
                span_data["UDL Length"] = load_entries["length"].get()

            ListofSpans[span_index] = span_data
            write_json("SpanData.json", ListofSpans)
            print(ListofSpans)

        SaveButton = ttk.Button(frame, text="Save", command=save_span)
        SaveButton.grid(row=5, column=0)

    SpanProceedButton = ttk.Button(beam_frame,text="Proceed",command=SpanProceed,)
    SpanProceedButton.grid(row=1, column=2)


# DRAWING HELPERS

def draw_support_symbol(canvas, x, y, support):
    if support == "Fixed Support":
        canvas.create_rectangle(x - 10, y - 45, x + 10, y + 45, fill="black")
    elif support == "Hinge Support":
        canvas.create_polygon(x, y, x - 20, y + 30, x + 20, y + 30, fill="gray")
    elif support == "Roller Support":
        canvas.create_polygon(x, y, x - 20, y + 30, x + 20, y + 30, fill="gray")
        canvas.create_oval(x - 15, y + 30, x - 5, y + 40)
        canvas.create_oval(x + 5, y + 30, x + 15, y + 40)
    else:
        canvas.create_line(x, y, x, y + 30, width=2)


def draw_span_loads(canvas, span, left_x, y, scale):
    load_type = span.get("Load Type", "No Load")

    if load_type == "Point Load":
        P = safe_float(span.get("Point Load Value"), 0.0)
        a = safe_float(span.get("Point Load Position"), 0.0)
        load_x = left_x + a * scale
        canvas.create_line(load_x, y - 80, load_x, y, arrow=tk.LAST, width=2)
        canvas.create_text(load_x, y - 95, text=f"{fmt(P)} kN")

    elif load_type == "UDL":
        w = safe_float(span.get("UDL Intensity"), 0.0)
        start = safe_float(span.get("UDL Start"), 0.0)
        length = safe_float(span.get("UDL Length"), span_length(span))
        x1 = left_x + start * scale
        x2 = left_x + (start + length) * scale
        arrow_count = max(3, int((x2 - x1) / 25))

        for i in range(arrow_count + 1):
            x = x1 + (x2 - x1) * i / arrow_count
            canvas.create_line(x, y - 60, x, y, arrow=tk.LAST)

        canvas.create_line(x1, y - 60, x2, y - 60, width=2)
        canvas.create_text((x1 + x2) / 2, y - 78, text=f"{fmt(w)} kN/m")


def draw_overhang_loads(canvas, side, support_x, y, scale):
    length = overhang_length(side)
    if length <= 0:
        return

    for load in normalized_overhang_loads(side):
        direction = -1 if side == "Left" else 1

        if load["kind"] == "point":
            x = support_x + direction * load["d"] * scale
            canvas.create_line(x, y - 75, x, y, arrow=tk.LAST, width=2)
            canvas.create_text(x, y - 90, text=f"{fmt(load['P'])} kN")

        elif load["kind"] == "udl":
            x1 = support_x + direction * load["start"] * scale
            x2 = support_x + direction * load["end"] * scale
            if x2 < x1:
                x1, x2 = x2, x1

            arrow_count = max(3, int(abs(x2 - x1) / 25))
            for i in range(arrow_count + 1):
                x = x1 + (x2 - x1) * i / arrow_count
                canvas.create_line(x, y - 55, x, y, arrow=tk.LAST)
            canvas.create_line(x1, y - 55, x2, y - 55, width=2)
            canvas.create_text((x1 + x2) / 2, y - 72, text=f"{fmt(load['w'])} kN/m")


def draw_beam():
    if app is None:
        return

    try:
        validate_model()
    except ValueError as exc:
        show_error("Beam Drawing", str(exc))
        return

    drawWindow = tk.Toplevel()
    drawWindow.title("Continuous Beam Drawing")
    drawWindow.geometry("1200x600")

    canvas = Canvas(drawWindow, bg="white")
    canvas.pack(fill="both", expand=True)

    geometry = beam_geometry()
    total_length = max(geometry["total_length"], 1.0)
    scale = min(90, 950 / total_length)
    start_x = 100
    y = 300

    support_positions = [start_x + position * scale for position in geometry["support_positions"]]

    beam_start = start_x
    beam_end = start_x + geometry["total_length"] * scale

    canvas.create_line(beam_start, y, beam_end, y, width=5)

    left_support_x = support_positions[0]
    if geometry["left_overhang"] > 0:
        canvas.create_text((beam_start + left_support_x) / 2, y + 35,text=f"Left overhang {fmt(geometry['left_overhang'])} m",)
        draw_overhang_loads(canvas, "Left", left_support_x, y, scale)

    current_x = left_support_x
    for i, span in enumerate(ListofSpans):
        L = span_length(span)
        end_x = current_x + L * scale
        canvas.create_text((current_x + end_x) / 2, y - 20, text=f"Span {i + 1}")
        canvas.create_text((current_x + end_x) / 2, y + 35, text=f"{fmt(L)} m")
        draw_span_loads(canvas, span, current_x, y, scale)
        current_x = end_x

    if geometry["right_overhang"] > 0:
        right_support_x = support_positions[-1]
        canvas.create_text(
            (right_support_x + beam_end) / 2,y + 35,text=f"Right overhang {fmt(geometry['right_overhang'])} m",)
        draw_overhang_loads(canvas, "Right", right_support_x, y, scale)

    for i, x in enumerate(support_positions):
        draw_support_symbol(canvas, x, y, support_type(i))
        canvas.create_text(x, y + 60, text=f"S{i + 1}")





def import_matplotlib():
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:
        raise RuntimeError( "Matplotlib is required for plotted SFD/BMD diagrams.\n\n"
            "Install it with:\n"
            "python -m pip install matplotlib"
        ) from exc

    return plt


def plot_diagram(kind):
    """Plot SFD/BMD with Matplotlib using solved three-moment results."""

    try:
        result = analyze_beam()
        plt = import_matplotlib()
    except (RuntimeError, ValueError) as exc:
        show_error(kind, str(exc))
        return

    geometry = beam_geometry()
    points = diagram_points(kind)
    x_values = [point[0] for point in points]
    y_values = [point[1] for point in points]
    support_positions = geometry["support_positions"]

    unit = "kN" if kind == "SFD" else "kNm"
    ylabel = "Shear Force (kN)" if kind == "SFD" else "Bending Moment (kNm)"
    title = ("Shear Force Diagram from Solved Reactions"
        if kind == "SFD"
        else "Final Bending Moment Diagram from Three-Moment Support Moments"
    )

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.axhline(0, color="black", linewidth=1.0)
    ax.plot(x_values, y_values, color="#0b5cad", linewidth=2.2)
    ax.fill_between(x_values, y_values, 0, color="#0b5cad", alpha=0.14)

    for i, position in enumerate(support_positions):
        ax.axvline(position, color="#777777", linestyle="--", linewidth=0.8, alpha=0.55)
        ax.text(position,0,f" S{i + 1}",rotation=90,va="bottom",ha="left",fontsize=9,color="#333333",)

    if kind == "SFD":
        for i, reaction in enumerate(result["reactions"]):
            ax.annotate(
                f"R{i + 1}={fmt(reaction, 5)} kN",
                xy=(support_positions[i], 0),
                xytext=(4, 18 + 12 * (i % 2)),
                textcoords="offset points",
                fontsize=9,
                arrowprops={"arrowstyle": "-", "color": "#555555", "lw": 0.8},)
    else:
        for i, moment in enumerate(result["moments"]):
            ax.scatter([support_positions[i]],[moment],color="#b3261e",zorder=4,)
            ax.annotate(f"M{i + 1}={fmt(moment, 5)} kNm",xy=(support_positions[i], moment),xytext=(5, 12 if moment >= 0 else -20),textcoords="offset points",fontsize=9,color="#b3261e",)

    ax.set_title(title)
    ax.set_xlabel("Distance along beam (m)")
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", linestyle=":", linewidth=0.8)
    ax.set_xlim(0, max(geometry["total_length"], 1.0))

    max_abs = max([abs(value) for value in y_values] + [1.0])
    ax.set_ylim(-1.18 * max_abs, 1.18 * max_abs)
    ax.text(0.99,0.02,f"Max abs = {fmt(max_abs, 6)} {unit}", transform=ax.transAxes, ha="right", va="bottom",fontsize=9, bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "#bbbbbb"},)

    fig.tight_layout()
    plt.show()


def plot_sfd_bmd():
    """Plot SFD and final BMD together after solving the three-moment equations."""

    try:
        result = analyze_beam()
        plt = import_matplotlib()
    except (RuntimeError, ValueError) as exc:
        show_error("SFD/BMD Plot", str(exc))
        return

    geometry = beam_geometry()
    support_positions = geometry["support_positions"]
    sfd_points = diagram_points("SFD")
    bmd_points = diagram_points("BMD")

    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    plot_specs = [(axes[0], sfd_points, "Shear Force Diagram", "Shear Force (kN)", "#0b5cad"),(axes[1],bmd_points,"Final Bending Moment Diagram","Bending Moment (kNm)","#b3261e",),
    ]

    for ax, points, title, ylabel, color in plot_specs:
        x_values = [point[0] for point in points]
        y_values = [point[1] for point in points]
        max_abs = max([abs(value) for value in y_values] + [1.0])

        ax.axhline(0, color="black", linewidth=1.0)
        ax.plot(x_values, y_values, color=color, linewidth=2.2)
        ax.fill_between(x_values, y_values, 0, color=color, alpha=0.14)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", linestyle=":", linewidth=0.8)
        ax.set_ylim(-1.18 * max_abs, 1.18 * max_abs)

        for i, position in enumerate(support_positions):
            ax.axvline(position, color="#777777",linestyle="--",linewidth=0.8,alpha=0.55)
            ax.text(position,0, f" S{i + 1}",rotation=90, va="bottom",ha="left",fontsize=9,color="#333333",)

    for i, reaction in enumerate(result["reactions"]):
        axes[0].annotate( f"R{i + 1}={fmt(reaction, 5)} kN",xy=(support_positions[i], 0), xytext=(4, 18 + 12 * (i % 2)), textcoords="offset points",fontsize=9,arrowprops={"arrowstyle": "-", "color": "#555555", "lw": 0.8},)

    for i, moment in enumerate(result["moments"]):
        axes[1].scatter([support_positions[i]], [moment], color="#b3261e", zorder=4)
        axes[1].annotate(f"M{i + 1}={fmt(moment, 5)} kNm",xy=(support_positions[i], moment), xytext=(5, 12 if moment >= 0 else -20),
            textcoords="offset points",fontsize=9,color="#b3261e",)

    axes[1].set_xlabel("Distance along beam (m)")
    axes[1].set_xlim(0, max(geometry["total_length"], 1.0))
    fig.suptitle("SFD and BMD from Three-Moment Analysis", fontsize=14)
    fig.tight_layout()
    plt.show()


def draw_sfd():
    plot_diagram("SFD")


def draw_bmd():
    plot_diagram("BMD")


#  RESULT WINDOWS 
def show_moments():
    if not SolvedMoments:
        show_info("Support Moments", "Run Three Moment Equation first.")
        return

    top = tk.Toplevel()
    top.title("Support Moments")
    top.geometry("450x500")

    frame = tk.LabelFrame(top, text="Calculated Support Moments")
    frame.pack(padx=20, pady=20, fill="both", expand=True)

    for i, value in enumerate(SolvedMoments):
        lbl = ttk.Label(frame, text=f"Support M{i + 1}")
        lbl.grid(row=i, column=0, padx=10, pady=10)

        ent = ttk.Entry(frame, width=25)
        ent.grid(row=i, column=1, padx=10, pady=10)
        ent.insert(0, f"{value} kNm")


def show_reactions():
    if not SolvedReactions:
        show_info("Support Reactions", "Run Three Moment Equation first.")
        return

    top = tk.Toplevel()
    top.title("Support Reactions")
    top.geometry("450x500")

    frame = tk.LabelFrame(top, text="Calculated Support Reactions")
    frame.pack(padx=20, pady=20, fill="both", expand=True)

    for i, value in enumerate(SolvedReactions):
        lbl = ttk.Label(frame, text=f"Support R{i + 1}")
        lbl.grid(row=i, column=0, padx=10, pady=10)

        ent = ttk.Entry(frame, width=25)
        ent.grid(row=i, column=1, padx=10, pady=10)
        ent.insert(0, f"{value} kN")


def show_equations():
    if not EquationRows:
        show_info("Three-Moment Equations", "Run Three Moment Equation first.")
        return

    if app is None:
        return

    top = tk.Toplevel()
    top.title("Symbolic Three-Moment Equations")
    top.geometry("900x650")

    text = tk.Text(top, wrap="word")
    text.pack(fill="both", expand=True, padx=10, pady=10)

    text.insert("end","Sign convention: sagging moments are positive; hogging moments are negative.\n\n",)

    text.insert("end", "THREE-MOMENT EQUATIONS\n\n")
    for row in EquationRows:
        text.insert("end", f"{row['name']}\n")
        text.insert("end", f"{row['symbol']}\n")
        text.insert("end", f"Numeric: {make_numeric_equation(row['coeffs'], row['rhs'])}\n\n")

    text.insert("end", "A*xbar TERMS\n\n")
    for i, stats in enumerate(AnalysisResult.get("span_stats", [])):
        text.insert("end",(f"Span {i + 1}: " f"A={fmt(stats['A'], 6)}, "f"A*xbar_left={fmt(stats['Ax_left'], 6)}, "f"A*xbar_right={fmt(stats['Ax_right'], 6)}\n"),)

    text.insert("end", "\nSOLVED SUPPORT MOMENTS\n\n")
    for i, value in enumerate(SolvedMoments):
        text.insert("end", f"M{i + 1} = {value} kNm\n")

    text.insert("end", "\nSOLVED SUPPORT REACTIONS\n\n")
    for i, value in enumerate(SolvedReactions):
        text.insert("end", f"R{i + 1} = {value} kN\n")


#  APP MENU 

def build_app():
    global app, menuBar

    app = tk.Tk()
    app.title("Beam Solver")
    app.geometry("1000x1000")
    app.configure(bg="lightblue")
    app.resizable(True, True)

    menuBar = tk.Menu(app)
    app.configure(menu=menuBar)

    fileMenu = tk.Menu(menuBar, tearoff=0)
    menuBar.add_cascade(label="File", menu=fileMenu)
    fileMenu.add_command(label="New file..", command=newModel)

    analysisMenu = tk.Menu(menuBar, tearoff=0)
    menuBar.add_cascade(label="Analysis", menu=analysisMenu)
    analysisMenu.add_command(label="End Moments / Overhangs",command=get_end_support_moments,)
    analysisMenu.add_command(label="Three Moment Equation",command=three_moment_equation,)
    analysisMenu.add_command(label="Show Moments", command=show_moments)
    analysisMenu.add_command(label="Show Reactions", command=show_reactions)
    analysisMenu.add_command(label="Show Equations", command=show_equations)

    DrawMenu = tk.Menu(menuBar, tearoff=0)
    menuBar.add_cascade(label="Draw", menu=DrawMenu)
    DrawMenu.add_command(label="Draw Continuous Beam", command=draw_beam)
    DrawMenu.add_command(label="Plot SFD", command=draw_sfd)
    DrawMenu.add_command(label="Plot Final BMD", command=draw_bmd)
    DrawMenu.add_command(label="Plot SFD + Final BMD", command=plot_sfd_bmd)

    newModel()
    return app


if __name__ == "__main__":
    build_app()
    app.mainloop()
