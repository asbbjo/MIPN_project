using UnityEngine;
using System.Collections.Generic;

public class RRT_drawlines : MonoBehaviour
{
    //public Text statusText; // UI text to display status (optional)
    public GameObject startObject; // Start waypoint assigned in the Inspector
    public GameObject endObject; // End waypoint assigned in the Inspector
    public GameObject boundaryObstacle; // Boundary assigned in the Inspector
    public List<GameObject> obstacles; // List of all obstacle objects assigned in the Inspector
    private List<GameObject> newPoints = new List<GameObject>();
    private Dictionary<GameObject, Dictionary<GameObject, float>> graph = new Dictionary<GameObject, Dictionary<GameObject, float>>();

    void Start()
    {
        Debug.Log("Start: Checking assignments before initialization");
        Debug.Log("startObject: " + (startObject != null ? startObject.name : "null"));
        Debug.Log("endObject: " + (endObject != null ? endObject.name : "null"));

        AddMeshCollidersToObstacles();
        AddMeshCollidersToStartEndObjects();

        if (startObject == null || endObject == null)
        {
            Debug.LogError("Start or End waypoint is not assigned. Please assign waypoints in the Inspector.");
            return;
        }

        GenerateRRTPoints(startObject, endObject, 100, 3, 6); //#of points, #distance and #number of branches
        Debug.Log("RRTPoints found."); 

        InitializeGraph();
        Debug.Log("Graph initialized.");


        // Automatically find the path at the start
        FindPath(startObject, endObject);
    }

    void AddMeshCollidersToObstacles()
    {
        foreach (var obstacle in obstacles)
        {
            MeshCollider meshCollider_i = obstacle.GetComponent<MeshCollider>();
            if (meshCollider_i == null)
            {
                meshCollider_i = obstacle.AddComponent<MeshCollider>();
                Debug.Log("Mesh Collider added to " + obstacle.name);
            }
        }

        MeshCollider meshCollider = boundaryObstacle.GetComponent<MeshCollider>();
        if (meshCollider == null)
        {
            meshCollider = boundaryObstacle.AddComponent<MeshCollider>();
            Debug.Log("Mesh Collider (for boundary) added to " + boundaryObstacle.name);
        }
    }

    void AddMeshCollidersToStartEndObjects()
    {
        List<GameObject> objects = new List<GameObject> { startObject, endObject };

        foreach (var my_object in objects)
        {
            if (my_object != null)
            {
                MeshCollider meshCollider = my_object.GetComponent<MeshCollider>();
                if (meshCollider == null)
                {
                    meshCollider = my_object.AddComponent<MeshCollider>();
                    Debug.Log("Mesh Collider added to " + my_object.name);
                }
            }
            else
            {
                Debug.LogWarning("Point is null.");
            }
        }
    }

    void AddMeshColliderToObject(GameObject new_object)
    {
        MeshCollider meshCollider = new_object.GetComponent<MeshCollider>();
        if (meshCollider == null)
        {
            meshCollider = new_object.AddComponent<MeshCollider>();
            Debug.Log("Mesh Collider added to " + new_object.name);
        }
    }

   bool IsObstacleBetween(GameObject object1, GameObject object2)
    {
        Vector3 direction = object2.transform.position - object1.transform.position;
        float distance = direction.magnitude;
        Ray ray = new Ray(object1.transform.position, direction);

        foreach (var obstacle in obstacles)
        {
            // Ensure the obstacle has a collider
            Collider collider = obstacle.GetComponent<Collider>();
            if (collider == null) continue;

            // Check if the ray intersects the collider
            if (collider.Raycast(ray, out RaycastHit hitInfo, distance))
            {
                // Check if the closest point on the collider is within the distance
                Vector3 closestPointToObstacle = collider.ClosestPoint(object1.transform.position);
                if (Vector3.Distance(object1.transform.position, closestPointToObstacle) < distance)
                {   
                    return true;
                }
            }
        }
        return false;
    }


    void GenerateRRTPoints(GameObject start, GameObject target, int numberOfPoints, float maxDistance, int numberOfBranches)
    {
        GameObject currentPoint = start;

        for (int i = 0; i < numberOfPoints; i++)
        {
            Vector3 newPoint = new Vector3(
                Random.Range(currentPoint.transform.position.x - maxDistance, currentPoint.transform.position.x + maxDistance),
                Random.Range(currentPoint.transform.position.y - maxDistance, currentPoint.transform.position.y + maxDistance),
                Random.Range(currentPoint.transform.position.z - maxDistance, currentPoint.transform.position.z + maxDistance)
            );

            GameObject newObject = new GameObject("New Waypoint" + i);
            newObject.transform.position = newPoint;


            for (int j = 0; j < numberOfBranches; j++)
            {
                Vector3 iterativePoint = new Vector3(
                    Random.Range(currentPoint.transform.position.x - maxDistance, currentPoint.transform.position.x + maxDistance),
                    Random.Range(currentPoint.transform.position.y - maxDistance, currentPoint.transform.position.y + maxDistance),
                    Random.Range(currentPoint.transform.position.z - maxDistance, currentPoint.transform.position.z + maxDistance)
                );

                GameObject iterativeObject = new GameObject("New Waypoint" + i);
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

            if (!IsObstacleBetween(currentPoint, newObject))
            {
                float distanceToCurrent = Vector3.Distance(currentPoint.transform.position, newObject.transform.position);

                if (!graph.ContainsKey(currentPoint))
                {
                    graph[currentPoint] = new Dictionary<GameObject, float>();
                }
                if (!graph.ContainsKey(newObject))
                {
                    graph[newObject] = new Dictionary<GameObject, float>();
                }

                graph[currentPoint][newObject] = distanceToCurrent;
                graph[newObject][currentPoint] = distanceToCurrent;

                currentPoint = newObject;
                newPoints.Add(newObject);
                Debug.Log("Added new waypoint " + newObject.name + " closer to target.");
            }
            else
            {
                Destroy(newObject);
                Debug.Log("New waypoint " + i + " is not valid due to obstacles.");
            }

            // Check if the current point is close enough to the target
            /*if (Vector3.Distance(currentPoint.transform.position, target.transform.position) < maxDistance)
            {
                float distanceToTarget = Vector3.Distance(currentPoint.transform.position, target.transform.position);
                graph[currentPoint][target] = distanceToTarget;
                graph[target][currentPoint] = distanceToTarget;
                Debug.Log("Connected to target waypoint.");
                break;
            }*/
        }
    }

    void InitializeGraph()
    {
        Debug.Log("Initializing graph...");
        graph[startObject] = new Dictionary<GameObject, float>();
        graph[endObject] = new Dictionary<GameObject, float>();

        foreach (var point in newPoints)
        {
            if (!graph.ContainsKey(point))
            {
                graph[point] = new Dictionary<GameObject, float>();
            }
        }

        // Connect startObject to the first new point
        if (newPoints.Count > 0 && !IsObstacleBetween(startObject, newPoints[0]))
        {
            float distance = Vector3.Distance(startObject.transform.position, newPoints[0].transform.position);
            graph[startObject][newPoints[0]] = distance;
            graph[newPoints[0]][startObject] = distance;
            Debug.Log("Added edge from " + startObject.name + " to " + newPoints[0].name + " with distance " + distance);
        }

        // Connect all new points to each other sequentially
        for (int i = 0; i < newPoints.Count - 1; i++)
        {
            if (!IsObstacleBetween(newPoints[i], newPoints[i + 1]))
            {
                float distance = Vector3.Distance(newPoints[i].transform.position, newPoints[i + 1].transform.position);
                graph[newPoints[i]][newPoints[i + 1]] = distance;
                graph[newPoints[i + 1]][newPoints[i]] = distance;
                Debug.Log("Added edge from " + newPoints[i].name + " to " + newPoints[i + 1].name + " with distance " + distance);
            }
        }

        // Connect the last new point to endObject
        if (newPoints.Count > 0 && !IsObstacleBetween(newPoints[newPoints.Count - 1], endObject))
        {
            float distance = Vector3.Distance(newPoints[newPoints.Count - 1].transform.position, endObject.transform.position);
            graph[newPoints[newPoints.Count - 1]][endObject] = distance;
            graph[endObject][newPoints[newPoints.Count - 1]] = distance;
            Debug.Log("Added edge from " + newPoints[newPoints.Count - 1].name + " to " + endObject.name + " with distance " + distance);
        }
    }

    void FindPath(GameObject start, GameObject end)
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
        Debug.Log("Path found: " + string.Join(" -> ", path.ConvertAll(node => node.name).ToArray()));
        VisualizePath(path);
    }

    void VisualizePath(List<GameObject> path)
    {
        for (int i = 0; i < path.Count - 1; i++)
        {
            if (!IsObstacleBetween(path[i], path[i + 1]))
            {
                Debug.DrawLine(path[i].transform.position, path[i + 1].transform.position, Color.red, 25f);
                Debug.Log("Drawing line from " + path[i].name + " to " + path[i + 1].name);
            }
            else
            {
                Debug.Log("Obstacle detected between " + path[i].name + " and " + path[i + 1].name + ". Path not visualized.");
                return; // Stop visualization if an obstacle is detected
            }
        }
    }


}