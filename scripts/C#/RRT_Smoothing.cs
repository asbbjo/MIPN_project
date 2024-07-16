using UnityEngine;
using System.Collections.Generic;
using System;
using System.Diagnostics; // Necesario para usar Stopwatch

public class RRT : MonoBehaviour
{
    //public Text statusText; // UI text to display status (optional)
    public GameObject startObject; // Start waypoint assigned in the Inspector
    public GameObject endObject; // End waypoint assigned in the Inspector
    public GameObject boundaryObstacle; // Boundary assigned in the Inspector
    public List<GameObject> obstacles; // List of all obstacle objects assigned in the Inspector
    public int numberOfPoints;
    public float maxDistance;
    public int numberOfBranches;
    public int threshold;
    private List<GameObject> newPoints = new List<GameObject>();
    private List<Vector3> newPointsSmooth = new List<Vector3>();
    private List<GameObject> tree = new List<GameObject>();
    private Dictionary<GameObject, Dictionary<GameObject, float>> graph = new Dictionary<GameObject, Dictionary<GameObject, float>>();
    private Dictionary<GameObject, Dictionary<GameObject, float>> graphSmooth = new Dictionary<GameObject, Dictionary<GameObject, float>>();
    List<GameObject> PruebaObjetos = new List<GameObject>();
    List<Vector3> NewPointsSmoothed = new List<Vector3>();
    List<Vector3> NewPointsSmoothedBIS = new List<Vector3>();
    List<Vector3> NewPointsSmoothedFinal = new List<Vector3>();

    void Start()
    {
        UnityEngine.Debug.Log("Start: Checking assignments before initialization");
        UnityEngine.Debug.Log("startObject: " + (startObject != null ? startObject.name : "null"));
        UnityEngine.Debug.Log("endObject: " + (endObject != null ? endObject.name : "null"));

        AddMeshCollidersToObstacles();
        AddMeshCollidersToStartEndObjects();

        if (startObject == null || endObject == null)
        {
            UnityEngine.Debug.LogError("Start or End waypoint is not assigned. Please assign waypoints in the Inspector.");
            return;
        }

        graph[startObject] = new Dictionary<GameObject, float>();
        graph[endObject] = new Dictionary<GameObject, float>();

        GenerateRRTPoints(startObject, endObject, numberOfPoints, maxDistance, numberOfBranches); //#of points, #distance and #number of branches
        VisualizePath(newPoints, Color.red);
        newPoints = PathPruning (newPoints);
        VisualizePath(newPoints, Color.green);
        //NewPointsSmoothed = SmoothPathRRT1(newPoints, endObject, startObject);
        NewPointsSmoothedFinal=SmoothPathRRT1(newPoints, endObject, startObject);
        
        for(int l=0 ; l < NewPointsSmoothedFinal.Count ; l++){
            UnityEngine.Debug.Log("nodo nro:"+l);
            GameObject newObject = new GameObject("New Waypoint Smooth" + l);
            newObject.transform.position=NewPointsSmoothedFinal[l];
            PruebaObjetos.Add(newObject);
    }
        VisualizePath(PruebaObjetos, Color.blue);


    }

List<GameObject> PathPruning (List<GameObject> newPoints ){
    int N = newPoints.Count;
    List<GameObject> newPointsPruned = new List<GameObject>();
    List<GameObject> newPointsPrunedOK = new List<GameObject>();
    int j=N-1;
    newPointsPruned.Add(newPoints[j]);
    while(j>1){
        for(int i=1 ; i<j ; i++){
            if(!IsObstacleBetween(newPoints[j], newPoints[i])){
                    newPointsPruned.Add(newPoints[i]);
                    j=i;
            }
        }

    }
        for (int k = newPointsPruned.Count-1; k >= 0; k--) {
        newPointsPrunedOK.Add(newPointsPruned[k]);
    }
     UnityEngine.Debug.Log("newPointsPrunedOK:");
    foreach (GameObject obj in newPointsPrunedOK) {
        UnityEngine.Debug.Log(obj);
    }

    UnityEngine.Debug.Log("newPointsPruned:");
    foreach (GameObject obj in newPointsPruned) {
        UnityEngine.Debug.Log(obj);
    }
    return newPointsPrunedOK;}

void AddMeshCollidersToObstacles(){
        foreach (var obstacle in obstacles)
        {
            MeshCollider meshCollider_i = obstacle.GetComponent<MeshCollider>();
            if (meshCollider_i == null)
            {
                meshCollider_i = obstacle.AddComponent<MeshCollider>();
                UnityEngine.Debug.Log("Mesh Collider added to " + obstacle.name);
            }
        }

        MeshCollider meshCollider = boundaryObstacle.GetComponent<MeshCollider>();
        if (meshCollider == null)
        {
            meshCollider = boundaryObstacle.AddComponent<MeshCollider>();
            UnityEngine.Debug.Log("Mesh Collider (for boundary) added to " + boundaryObstacle.name);
        }}

void AddMeshCollidersToStartEndObjects(){
        List<GameObject> objects = new List<GameObject> { startObject, endObject };

        foreach (var my_object in objects)
        {
            if (my_object != null)
            {
                MeshCollider meshCollider = my_object.GetComponent<MeshCollider>();
                if (meshCollider == null)
                {
                    meshCollider = my_object.AddComponent<MeshCollider>();
                    UnityEngine.Debug.Log("Mesh Collider added to " + my_object.name);
                }
            }
            else
            {
                UnityEngine.Debug.LogWarning("Point is null.");
            }
        }}

void AddMeshColliderToObject(GameObject new_object){
        MeshCollider meshCollider = new_object.GetComponent<MeshCollider>();
        if (meshCollider == null)
        {
            meshCollider = new_object.AddComponent<MeshCollider>();
            UnityEngine.Debug.Log("Mesh Collider added to " + new_object.name);
        }}

bool IsObstacleBetween(GameObject object1, GameObject object2){
    Vector3 direction = (object2.transform.position - object1.transform.position);
    float distance = direction.magnitude;
    float threshold_obstacles = distance * threshold;
    Ray ray = new Ray(object1.transform.position, direction);

    foreach (var obstacle in obstacles){
        // Ensure the obstacle has a collider
        Collider collider = obstacle.GetComponent<Collider>();
        if (collider == null) continue;

        // Check if the ray intersects the collider
        if (collider.Raycast(ray, out RaycastHit hitInfo, threshold_obstacles)){
            // Check if the closest point on the collider is within the distance
            Vector3 closestPointToObstacle = collider.ClosestPoint(object1.transform.position);
            if (Vector3.Distance(object1.transform.position, closestPointToObstacle) < threshold_obstacles){   
                return true;
                }
        }
    }
    return false;}


void GenerateRRTPoints(GameObject start, GameObject target, int numberOfPoints, float maxDistance, int numberOfBranches){
    GameObject currentObject = start;
    for (int i = 0; i < numberOfPoints; i++){
        Vector3 newPoint = new Vector3(
        UnityEngine.Random.Range(currentObject.transform.position.x - maxDistance, currentObject.transform.position.x + maxDistance),
        UnityEngine.Random.Range(currentObject.transform.position.y - maxDistance, currentObject.transform.position.y + maxDistance),
        UnityEngine.Random.Range(currentObject.transform.position.z - maxDistance, currentObject.transform.position.z + maxDistance));

        GameObject newObject = new GameObject("New Waypoint " + i);
        newObject.transform.position = newPoint;

        for (int j = 0; j < numberOfBranches; j++){
            Vector3 iterativePoint = new Vector3(
            UnityEngine.Random.Range(currentObject.transform.position.x - maxDistance, currentObject.transform.position.x + maxDistance),
            UnityEngine.Random.Range(currentObject.transform.position.y - maxDistance, currentObject.transform.position.y + maxDistance),
            UnityEngine.Random.Range(currentObject.transform.position.z - maxDistance, currentObject.transform.position.z + maxDistance)
                );

                GameObject iterativeObject = new GameObject("New Waypoint " + i + "-" + j);
                iterativeObject.transform.position = iterativePoint;

                float newDist = Vector3.Distance(iterativePoint, target.transform.position);
                float oldDist = Vector3.Distance(newPoint, target.transform.position);

                if (newDist < oldDist)
                {
                    newObject = iterativeObject;
                    newPoint = iterativePoint;
                }
        }
            
        AddMeshColliderToObject(newObject);

        if (!IsObstacleBetween(currentObject, newObject)){
            float distanceToCurrent = Vector3.Distance(currentObject.transform.position, newObject.transform.position);

            if (!graph.ContainsKey(currentObject)){graph[currentObject] = new Dictionary<GameObject, float>();}
            if (!graph.ContainsKey(newObject)){graph[newObject] = new Dictionary<GameObject, float>();}

            graph[currentObject][newObject] = distanceToCurrent;
            graph[newObject][currentObject] = distanceToCurrent;

            currentObject = newObject;
            newPoints.Add(newObject);
        }
        else{
            Destroy(newObject);
            UnityEngine.Debug.Log("New waypoint " + i + "-X is not valid due to obstacles.");
        }

            // Check if the current point is close enough to the target
        if (Vector3.Distance(currentObject.transform.position, target.transform.position) < maxDistance/2){
            float distanceToTarget = Vector3.Distance(currentObject.transform.position, target.transform.position);
            graph[currentObject][target] = distanceToTarget;
            graph[target][currentObject] = distanceToTarget;
            //UnityEngine.Debug.Log("Connected to target waypoint.");
            break;
        }
    }
    newPoints.Add(target);}






    List<GameObject> FindPath(GameObject start, GameObject end)
    {
        Dictionary<GameObject, float> distances = new Dictionary<GameObject, float>();
        Dictionary<GameObject, GameObject> previous = new Dictionary<GameObject, GameObject>();
        List<GameObject> unvisited = new List<GameObject>(graph.Keys);

        foreach (var node in graph.Keys)
        {
            distances[node] = float.PositiveInfinity;
            previous[node] = null;
        }

        distances[start] = 0;

        while (unvisited.Count > 0)
        {
            GameObject currentNode = null;
            float smallestDistance = float.PositiveInfinity;

            foreach (var node in unvisited)
            {
                if (distances[node] < smallestDistance)
                {
                    smallestDistance = distances[node];
                    currentNode = node;
                }
            }

            if (currentNode == null)
                break;

            unvisited.Remove(currentNode);

            foreach (var neighbor in graph[currentNode].Keys)
            {
                float newDistance = distances[currentNode] + graph[currentNode][neighbor];

                if (newDistance < distances[neighbor])
                {
                    distances[neighbor] = newDistance;
                    previous[neighbor] = currentNode;
                }
            }
        }

        List<GameObject> path = new List<GameObject>();
        GameObject current = end;

        while (current != null)
        {
            path.Add(current);
            current = previous[current];
        }

        path.Reverse();
        UnityEngine.Debug.Log("Path found: " + string.Join(" -> ", path.ConvertAll(node => node.name).ToArray()));

        return path;
    }


void VisualizePath(List<GameObject> path, Color color){
    for (int i = 0; i < path.Count - 1; i++){
        UnityEngine.Debug.DrawLine(path[i].transform.position, path[i + 1].transform.position, color, 5000f);
        UnityEngine.Debug.Log("Drawing line from " + path[i].name + " to " + path[i + 1].name);
    }
}


Vector3 CrossProduct(Vector3 v1, Vector3 v2){
        float X = v1.y * v2.z - v1.z * v2.y;
        float Y = v1.z * v2.x - v1.x * v2.z;
        float Z = v1.x * v2.y - v1.y * v2.x;

        return new Vector3(X,Y,Z);}


List<Vector3> SmoothPathRRT1(List<GameObject> waypoints, GameObject endObject, GameObject startObject){
    List<Vector3> smoothedPath = new List<Vector3>();
    List<Vector2> smoothedPath2D = new List<Vector2>();
    List<Vector3> NewPointsSmoothed = new List<Vector3>();
    List<GameObject> NewObjectsSmoothed = new List<GameObject>();
    GameObject ObjectA = new GameObject("ObjectA");
    GameObject ObjectB = new GameObject("ObjectB");
    List<Vector3> Prueba3D = new List<Vector3>();
    List<float> z3D = new List<float>();
    List<float> z3DPOINTS = new List<float>();
    List<GameObject> PruebaObjetos = new List<GameObject>();
    bool HitsObstacle = false;
    int count=0;
    List<Vector2> Points2D = new List<Vector2>();
    List<Vector3> Points3D = new List<Vector3>();

    Vector3 P0_3D = new Vector3();
    Vector3 P2_3D = new Vector3();
    Vector3 P3_3D = new Vector3();

    Vector3 ux = new Vector3();
    Vector3 uy = new Vector3();
    Vector3 uz = new Vector3();

    Vector4 P0_4D = new Vector4();
    Vector4 P2_4D = new Vector4();
    Vector4 P3_4D = new Vector4();
        
    P0_3D=waypoints[0].transform.position;
    for (int i = 1; i <= waypoints.Count - 2; i=i+1){
        UnityEngine.Debug.Log("iteracion " + i);
        P2_3D=waypoints[i].transform.position;
        P3_3D=waypoints[i+1].transform.position;

        ux = (P2_3D - P0_3D).normalized;
        uy = (P3_3D - P2_3D).normalized;
        uz = CrossProduct(ux, uy);
        uy = CrossProduct(ux,uz);

        float[,] TM = {
            { ux.x , uy.x , uz.x , P0_3D.x },
            { ux.y , uy.y , uz.y , P0_3D.y },
            { ux.z , uy.z , uz.z , P0_3D.z },
            {   0 ,    0 ,    0 ,    1     }
        };

        float[,] TM_inv = CalcularInversa(TM);
        P0_4D = new Vector4 (P0_3D.x, P0_3D.y, P0_3D.z, 1f);
        P2_4D = new Vector4(P2_3D.x, P2_3D.y, P2_3D.z, 1f);
        P3_4D = new Vector4(P3_3D.x, P3_3D.y, P3_3D.z, 1f);
        Vector4 MP0_4D = MultMatrix_vector(TM_inv , P0_4D);
        Vector4 MP2_4D = MultMatrix_vector(TM_inv , P2_4D);
        Vector4 MP3_4D = MultMatrix_vector(TM_inv , P3_4D);

        UnityEngine.Debug.Log("Puntos MP0_3D: " + MP0_4D.x + "---" + MP0_4D.y + "---"+MP0_4D.z + "---"+MP0_4D.w);
        UnityEngine.Debug.Log("Puntos MP2_3D: " + MP2_4D.x + "---" + MP2_4D.y + "---"+MP2_4D.z + "---"+MP2_4D.w);
        UnityEngine.Debug.Log("Puntos MP3_3D: " + MP3_4D.x + "---" + MP3_4D.y + "---"+MP3_4D.z + "---"+MP3_4D.w);

        
        Vector2 P0 =new Vector2 (MP0_4D.x, MP0_4D.y);
        Vector2 P2 = new Vector2(MP2_4D.x, MP2_4D.y);
        Vector2 P3 = new Vector2(MP3_4D.x, MP3_4D.y);

        float titadeg =180f - Vector2.Angle((P0-P2),(P3-P2));
        float tita = titadeg * Mathf.PI /180;
        float beta = tita/2;
        float k = Vector3.Distance(P2,P3);
        float h = k * (4.58f) / (6 * Mathf.Cos(tita));
        float g = 0.58f * h;
        
        float kMax = 2 * h * Mathf.Sin(tita)/(3*k*k);
        float dMax = 1.1228f * Mathf.Sin(beta)/(kMax*Mathf.Cos(beta)*Mathf.Cos(beta));
     

        Vector2 u1 = (P0 - P2).normalized; //P2P0 
        Vector2 u2 = (P3 - P2).normalized; //P2P3
        Vector2 P1 = P0 - u1 * g;
        float D32=(P3-P2).magnitude;

        float D02=(P0-P2).magnitude;
        UnityEngine.Debug.Log("medida D32= " + D32);
        UnityEngine.Debug.Log("medida D02= " + D02);
        float d = 0f;
        if(D32<D02){
            d=D32/2;
        }
        else{
            d=D02/2;
        }
    UnityEngine.Debug.Log("medida d= " + d);
    //float d = 0.03f*D.magnitude/Mathf.Sin(tita);
    count=0;
    do{
    HitsObstacle = false;
    /*float hb = 0.346f*d;
    float gb = 0.58f*hb;
    float kb = 1.31f * hb * Mathf.Cos(beta);
    float he = hb;
    float ge = 0.58f * hb;
    float ke = 1.31f * hb * Mathf.Cos(beta);
*/
    //u1 = (P2-P1).normalized;

        Vector2 B0 = P0; //P2 + d * u1; //P0
        Vector2 B1 = P2 + (D02/2) * u1;
        Vector2 B2 = B1 - (D02/4) * u1;
        
        Vector2 E3 = P3- (D32/2) * u2;
        Vector2 E2 = P2 + (D32/4) * u2;
        Vector2 E1 = E2 - (D32/8) * u2;

        Vector2 ud = (E1-B2).normalized;
        Vector2 B3 = B2 + ud * (ud.magnitude)/2;

        Vector2 E0 = E1 - ud * (ud.magnitude)/2;

            // Generate points along the first Bezier curve
            List<Vector2> curve1Points = GenerateBezierCurve(B0, B1, B2,B3); //new List<Vector2>{B0, B1, B2, B3};//
            List<Vector2> BPoints = new List<Vector2>{B0, B1, B2, B3};//
            // Generate points along the second Bezier curve
            List<Vector2> curve2Points =GenerateBezierCurve(E0, E1, E2, E3);// new List<Vector2>{E0, E1, E2, E3}; //
            List<Vector2> EPoints = new List<Vector2>{E0, E1, E2, E3}; //
            // Add the points of curve1Points to smoothedPath (excluding duplicates)
            foreach (var point in curve1Points){
                if (!smoothedPath2D.Contains(point)){smoothedPath2D.Add(point);}
            }

            // Add the points of curve2Points to smoothedPath (excluding duplicates)
            foreach (var point in curve2Points){
                if (!smoothedPath2D.Contains(point)){smoothedPath2D.Add(point);}
            }

            foreach (var point in BPoints){
                if (!Points2D.Contains(point)){Points2D.Add(point);}
            }

            // Add the points of curve2Points to smoothedPath (excluding duplicates)
            foreach (var point in EPoints){
                if (!Points2D.Contains(point)){Points2D.Add(point);}
            }

            smoothedPath =  ConvertTo3D(smoothedPath2D, TM, "xy");
            Points3D=ConvertTo3D(Points2D, TM, "xy");
            smoothedPath2D.Clear();
            Points2D.Clear();

            for(int j = 0; j < smoothedPath.Count-1; j++){
                ObjectA.transform.position=smoothedPath[j];
                ObjectB.transform.position=smoothedPath[j + 1];
                UnityEngine.Debug.Log("Distancia entre puntos: "+ (Vector3.Distance(smoothedPath[j],smoothedPath[j+1])));
                if((Vector3.Distance(smoothedPath[j],smoothedPath[j+1]))>10){
                    if(IsObstacleBetween(ObjectA, ObjectB)){
                        UnityEngine.Debug.Log("Entra en obstacleBetween para iteracion " + i);
                        d = d*0.2f; 
                        HitsObstacle=true;
                        count++;
                        if(count == 16){
                           //d=0.03f*d;
                        }
                    }
                }
                else{UnityEngine.Debug.Log("Distancia entre puntos: "+ (Vector3.Distance(smoothedPath[j],smoothedPath[j+1])));}
                   
            }
        }while(HitsObstacle && count <16);

        if(count <16){
            UnityEngine.Debug.Log("Se usa d= "+ d+ "para iteracion " + i);
            NewPointsSmoothed.AddRange(smoothedPath);
            P0_3D=NewPointsSmoothed[NewPointsSmoothed.Count-1];
            UnityEngine.Debug.Log("cantidad de nodos en NewPointsSmoothed:"+ NewPointsSmoothed.Count);
        }
        else{
            UnityEngine.Debug.Log("No se encontro d para iteracion " + i);
            NewPointsSmoothed.AddRange(smoothedPath);
            //NewPointsSmoothed.AddRange(Points3D);
            P0_3D=NewPointsSmoothed[NewPointsSmoothed.Count-1];
        }
    }
    //PruebaObjetos.Add(endObject);
    return NewPointsSmoothed;}


Vector2 BezierCurvePoint(float t, Vector2 p0, Vector2 p1, Vector2 p2, Vector2 p3)
    {
    float u = 1 - t;
    float tt = t * t;
    float uu = u * u;
    float uuu = uu * u;
    float ttt = tt * t;

    Vector2 p = uuu * p0; // (1-t)^3 * P0
    p += 3 * uu * t * p1; // 3 * (1-t)^2 * t * P1
    p += 3 * u * tt * p2; // 3 * (1-t) * t^2 * P2
    p += ttt * p3;        // t^3 * P3

    return p;
    }


Vector4 MultMatrix_vector(float[,] M,Vector4 P){
    Vector4 PN = new Vector4();
        PN.x = M[0, 0] * P.x + M[0, 1] * P.y + M[0, 2] * P.z+ M[0, 3] * P.w;
        PN.y = M[1, 0] * P.x + M[1, 1] * P.y +M[1, 2] * P.z+ M[1, 3] * P.w;
        PN.z = M[2, 0] * P.x + M[2, 1] * P.y + M[2, 2] * P.z+ M[2, 3] * P.w;
        PN.w = M[3, 0] * P.x + M[3, 1] * P.y + M[3, 2] * P.z+ M[3, 3] * P.w;
        return PN;
}

    // Function to generate points along a cubic Bezier curve
    List<Vector2> GenerateBezierCurve(Vector2 P0, Vector2 P1, Vector2 P2, Vector2 P3)
    {
        List<Vector2> curvePoints = new List<Vector2>();
        float tStep = 0.05f; // Step size for interpolation

        for (float t = 0; t <= 1; t += tStep)
        {
            Vector2 point = BezierCurvePoint(t,P0, P1, P2, P3);
            curvePoints.Add(point);
        }

        // Add the final point (P3)
        curvePoints.Add(P3);

        return curvePoints;
    }




List<Vector3> ConvertTo3D(List<Vector2> points2D,float[,] M, string plane)
{
    List<Vector3> points3DFinal = new List<Vector3>();
    Vector3 point3D = new Vector3();
    Vector4 point4D = new Vector4 ();
    int i=0;
    foreach (var point in points2D)
    {   UnityEngine.Debug.Log("Puntos 2D de ConvertTo3D: " + point3D.x + "," + point3D.y);
        if(plane=="xy"){point4D = new Vector4(point.x, point.y,0f,1f);}
        if(plane=="xz"){point4D = new Vector4(point.x,0f, point.y,1f);}
        if(plane=="yz"){point4D = new Vector4(0f,point.x,point.y,1f);}
        point4D = MultMatrix_vector(M, point4D);
        point3D= new Vector3(point4D.x, point4D.y,point4D.z);
        points3DFinal.Add(point3D);
        UnityEngine.Debug.Log("Puntos de ConvertTo3D: " + point3D.x + "," + point3D.y + "," + point3D.z);
        i++;
    }
    

    return points3DFinal;
}


static float[,] CalcularInversa(float[,] matriz)
    {
        int n = matriz.GetLength(0);
        float[,] adjunta = new float[n, n];
        float det = Determinante(matriz);

        if (det == 0)
        {
            throw new InvalidOperationException("La matriz no tiene inversa (es singular).");
        }

        float[,] cofactor = new float[n, n];

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                cofactor[i, j] = Mathf.Pow(-1, i + j) * Determinante(Submatriz(matriz, i, j));
            }
        }

        // Transponer la matriz de cofactores para obtener la adjunta
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                adjunta[j, i] = cofactor[i, j];
            }
        }

        // Calcular la matriz inversa
        float[,] inversa = new float[n, n];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                inversa[i, j] = adjunta[i, j] / det;
            }
        }

        return inversa;
    }

    static float Determinante(float[,] matriz)
    {
        int n = matriz.GetLength(0);
        if (n != matriz.GetLength(1))
        {
            throw new ArgumentException("La matriz debe ser cuadrada para calcular el determinante.");
        }

        if (n == 1)
        {
            return matriz[0, 0];
        }

        if (n == 2)
        {
            return matriz[0, 0] * matriz[1, 1] - matriz[0, 1] * matriz[1, 0];
        }

        float determinante = 0.0f;
        for (int j = 0; j < n; j++)
        {
            determinante += Mathf.Pow(-1, j) * matriz[0, j] * Determinante(Submatriz(matriz, 0, j));
        }

        return determinante;
    }

    static float[,] Submatriz(float[,] matriz, int filaEliminar, int columnaEliminar)
    {
        int n = matriz.GetLength(0);
        float[,] submatriz = new float[n - 1, n - 1];
        int filaDestino = 0;
        int columnaDestino = 0;

        for (int fila = 0; fila < n; fila++)
        {
            if (fila == filaEliminar)
                continue;

            columnaDestino = 0;

            for (int columna = 0; columna < n; columna++)
            {
                if (columna == columnaEliminar)
                    continue;

                submatriz[filaDestino, columnaDestino] = matriz[fila, columna];
                columnaDestino++;
            }

            filaDestino++;
        }

        return submatriz;
    }

}